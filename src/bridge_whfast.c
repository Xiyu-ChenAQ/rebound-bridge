/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Preserve WHFast's Jacobi state across synchronized Bridge stages. */
#include "bridge_internal.h"

static struct reb_bridge_whfast_cache* component_cache(
    struct reb_bridge* bridge, const struct reb_simulation* sim
) {
    if (sim == bridge->main_sim) return &bridge->main_whfast_cache;
    for (int i = 0; i < bridge->n_subsystems; i++) {
        if (sim == bridge->subsystems[i].sub_sim) return &bridge->subsystems[i].whfast_cache;
    }
    return NULL;
}

static int supports_jacobi_kicks(const struct reb_simulation* sim) {
    if (sim->integrator != REB_INTEGRATOR_WHFAST || sim->N <= 0 || sim->N_var != 0 ||
        (sim->N_active != -1 && sim->N_active != sim->N) ||
        sim->ri_whfast.coordinates != REB_WHFAST_COORDINATES_JACOBI ||
        sim->ri_whfast.kernel != REB_WHFAST_KERNEL_DEFAULT ||
        sim->ri_whfast.corrector != 0 || sim->ri_whfast.corrector2 != 0 ||
        sim->ri_whfast.keep_unsynchronized != 0 ||
        sim->additional_forces || sim->pre_timestep_modifications ||
        sim->post_timestep_modifications || sim->heartbeat ||
        sim->collision != REB_COLLISION_NONE || sim->boundary != REB_BOUNDARY_NONE) return 0;
    if (!(sim->particles[0].m > 0.0)) return 0;
    double mass = 0.0;
    for (int i = 0; i < sim->N; i++) {
        if (!isfinite(sim->particles[i].m) || sim->particles[i].m < 0.0) return 0;
        mass += sim->particles[i].m;
    }
    return isfinite(mass);
}

static int cache_matches(const struct reb_bridge_whfast_cache* cache, const struct reb_simulation* sim) {
    if (!cache || !cache->valid || cache->count != sim->N ||
        cache->jacobi != sim->ri_whfast.p_jh || !sim->ri_whfast.p_jh ||
        sim->ri_whfast.N_allocated != (unsigned int)sim->N ||
        !sim->ri_whfast.is_synchronized || sim->ri_whfast.recalculate_coordinates_this_timestep) return 0;
    for (int i = 0; i < sim->N; i++) {
        const struct reb_particle* p = &sim->particles[i];
        const struct reb_bridge_cartesian_state* s = &cache->state[i];
        if (p->m != s->m || p->x != s->x || p->y != s->y || p->z != s->z ||
            p->vx != s->vx || p->vy != s->vy || p->vz != s->vz) return 0;
    }
    return 1;
}

static void remember_state(struct reb_bridge_whfast_cache* cache, const struct reb_simulation* sim) {
    if (!cache) return;
    cache->valid = 0;
    if (!supports_jacobi_kicks(sim) || !sim->ri_whfast.is_synchronized ||
        sim->ri_whfast.recalculate_coordinates_this_timestep || !sim->ri_whfast.p_jh) return;
    if (cache->capacity < sim->N) {
        struct reb_bridge_cartesian_state* state = realloc(cache->state, sizeof(*state) * (size_t)sim->N);
        /* Missing a snapshot only disables reuse; it must not lose a kick. */
        if (!state) return;
        cache->state = state;
        cache->capacity = sim->N;
    }
    for (int i = 0; i < sim->N; i++) {
        const struct reb_particle* p = &sim->particles[i];
        struct reb_bridge_cartesian_state* s = &cache->state[i];
        s->m=p->m; s->x=p->x; s->y=p->y; s->z=p->z;
        s->vx=p->vx; s->vy=p->vy; s->vz=p->vz;
    }
    cache->count = sim->N;
    cache->jacobi = sim->ri_whfast.p_jh;
    cache->valid = 1;
}

static void synchronize_for_kick(struct reb_simulation* sim) {
    if (sim->integrator != REB_INTEGRATOR_WHFAST) {
        reb_simulation_synchronize(sim);
        return;
    }
    /* Commit synchronization even if the user normally keeps a staggered state.
     * The final half drift belongs to the last executed step, not nominal dt. */
    const unsigned int keep = sim->ri_whfast.keep_unsynchronized;
    const double dt = sim->dt;
    sim->ri_whfast.keep_unsynchronized = 0;
    if (!sim->ri_whfast.is_synchronized && sim->dt_last_done != 0.0) sim->dt = sim->dt_last_done;
    reb_simulation_synchronize(sim);
    sim->dt = dt;
    sim->ri_whfast.keep_unsynchronized = keep;
}

void reb_bridge_internal_prepare_cross_kick(struct reb_bridge* bridge) {
    synchronize_for_kick(bridge->main_sim);
    for (int i = 0; i < bridge->n_subsystems; i++) synchronize_for_kick(bridge->subsystems[i].sub_sim);
}

enum REB_STATUS reb_bridge_internal_integrate_component(
    struct reb_bridge* bridge, struct reb_simulation* sim, double target
) {
    struct reb_bridge_whfast_cache* cache = component_cache(bridge, sim);
    /* With no cross interaction, retain the component integrator's own call
     * semantics, including its safe-mode arithmetic and synchronization. */
    const int coupled = bridge->n_subsystems > 0 && bridge->main_sim->N > 1;
    const int reuse = coupled && bridge->whfast_jacobi_updates && supports_jacobi_kicks(sim);
    const unsigned int safe = sim->ri_whfast.safe_mode;
    const unsigned int keep = sim->ri_whfast.keep_unsynchronized;
    if (reuse) {
        if (!cache_matches(cache, sim)) sim->ri_whfast.recalculate_coordinates_this_timestep = 1;
        sim->ri_whfast.safe_mode = 0;
    } else {
        if (coupled && sim->integrator == REB_INTEGRATOR_WHFAST) {
            /* Without a reusable snapshot, Cartesian state is authoritative. */
            sim->ri_whfast.recalculate_coordinates_this_timestep = 1;
        }
        if (cache) cache->valid = 0;
    }
    /* Bridge stages expose a synchronized physical state at their boundaries. */
    if (sim->integrator == REB_INTEGRATOR_WHFAST && !safe) sim->ri_whfast.keep_unsynchronized = 0;
    const enum REB_STATUS status = reb_simulation_integrate(sim, target);
    sim->ri_whfast.safe_mode = safe;
    sim->ri_whfast.keep_unsynchronized = keep;
    if (status == REB_STATUS_SUCCESS && reuse) remember_state(cache, sim);
    else if (cache) cache->valid = 0;
    return status;
}

void reb_bridge_internal_apply_velocity_kick(
    struct reb_bridge* bridge, struct reb_simulation* sim,
    const struct reb_bridge_vec3* acc, double dt
) {
    if (dt == 0.0) return;
    struct reb_bridge_whfast_cache* cache = component_cache(bridge, sim);
    if (bridge->whfast_jacobi_updates && supports_jacobi_kicks(sim) && cache_matches(cache, sim)) {
        /* For fixed masses the velocity transformation is linear. Transform
         * only dv, leaving the already integrated Jacobi positions untouched. */
        struct reb_particle* jacobi = sim->ri_whfast.p_jh;
        double mass = sim->particles[0].m;
        struct reb_bridge_vec3 weighted_acc = vec_scale(acc[0], mass);
        for (int i = 1; i < sim->N; i++) {
            const struct reb_bridge_vec3 relative_acc = vec_sub(acc[i], vec_scale(weighted_acc, 1.0 / mass));
            jacobi[i].vx += dt * relative_acc.x;
            jacobi[i].vy += dt * relative_acc.y;
            jacobi[i].vz += dt * relative_acc.z;
            weighted_acc = vec_add(weighted_acc, vec_scale(acc[i], sim->particles[i].m));
            mass += sim->particles[i].m;
        }
        const struct reb_bridge_vec3 com_acc = vec_scale(weighted_acc, 1.0 / mass);
        /* 0 号存整体质心 */
        jacobi[0].vx += dt * com_acc.x;
        jacobi[0].vy += dt * com_acc.y;
        jacobi[0].vz += dt * com_acc.z;
        reb_particles_transform_jacobi_to_inertial_posvel(sim->particles, jacobi, sim->particles, sim->N, sim->N);
        remember_state(cache, sim);
    } else {
        for (int i = 0; i < sim->N; i++) {
            sim->particles[i].vx += dt * acc[i].x;
            sim->particles[i].vy += dt * acc[i].y;
            sim->particles[i].vz += dt * acc[i].z;
        }
        if (sim->integrator == REB_INTEGRATOR_WHFAST) sim->ri_whfast.recalculate_coordinates_this_timestep = 1;
        if (cache) cache->valid = 0;
    }
}
