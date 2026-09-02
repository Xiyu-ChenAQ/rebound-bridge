/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "bridge_internal.h"

/* Cross-kick workspace and plain Bridge coupling. */

int reb_bridge_internal_ensure_kick_workspace(struct reb_bridge* bridge, int sub_particles, int main_particles) {
    if (sub_particles > bridge->cap_sub_particles) {
        struct reb_bridge_vec3* internal_acc =
            (struct reb_bridge_vec3*)realloc(
                bridge->internal_acc,
                sizeof(*internal_acc) * (size_t)sub_particles
            );
        if (!internal_acc) return bridge_error("out of memory");
        bridge->internal_acc = internal_acc;

        struct reb_bridge_vec3* source_acc =
            (struct reb_bridge_vec3*)realloc(
                bridge->source_acc,
                sizeof(*source_acc) * (size_t)sub_particles
            );
        if (!source_acc) return bridge_error("out of memory");
        bridge->source_acc = source_acc;

        struct reb_bridge_vec3* fg_internal_acc =
            (struct reb_bridge_vec3*)realloc(
                bridge->fg_internal_acc,
                sizeof(*fg_internal_acc) * (size_t)sub_particles
            );
        if (!fg_internal_acc) return bridge_error("out of memory");
        bridge->fg_internal_acc = fg_internal_acc;

        struct reb_bridge_vec3* fd_internal_acc =
            (struct reb_bridge_vec3*)realloc(
                bridge->fd_internal_acc,
                sizeof(*fd_internal_acc) * (size_t)sub_particles
            );
        if (!fd_internal_acc) return bridge_error("out of memory");
        bridge->fd_internal_acc = fd_internal_acc;

        struct reb_bridge_vec3* saved_sub_pos =
            (struct reb_bridge_vec3*)realloc(
                bridge->saved_sub_pos,
                sizeof(*saved_sub_pos) * (size_t)sub_particles
            );
        if (!saved_sub_pos) return bridge_error("out of memory");
        bridge->saved_sub_pos = saved_sub_pos;
        bridge->cap_sub_particles = sub_particles;
    }

    if (main_particles > bridge->cap_main_particles) {
        struct reb_bridge_vec3* main_acc =
            (struct reb_bridge_vec3*)realloc(
                bridge->main_acc,
                sizeof(*main_acc) * (size_t)main_particles
            );
        if (!main_acc) return bridge_error("out of memory");
        bridge->main_acc = main_acc;

        struct reb_bridge_vec3* fg_main_acc =
            (struct reb_bridge_vec3*)realloc(
                bridge->fg_main_acc,
                sizeof(*fg_main_acc) * (size_t)main_particles
            );
        if (!fg_main_acc) return bridge_error("out of memory");
        bridge->fg_main_acc = fg_main_acc;

        struct reb_bridge_vec3* fd_main_acc =
            (struct reb_bridge_vec3*)realloc(
                bridge->fd_main_acc,
                sizeof(*fd_main_acc) * (size_t)main_particles
            );
        if (!fd_main_acc) return bridge_error("out of memory");
        bridge->fd_main_acc = fd_main_acc;

        struct reb_bridge_vec3* saved_main_pos =
            (struct reb_bridge_vec3*)realloc(
                bridge->saved_main_pos,
                sizeof(*saved_main_pos) * (size_t)main_particles
            );
        if (!saved_main_pos) return bridge_error("out of memory");
        bridge->saved_main_pos = saved_main_pos;
        bridge->cap_main_particles = main_particles;
    }

    return 0;
}

int reb_bridge_internal_check_subsystem_time_sync(
    struct reb_bridge* bridge,
    struct reb_bridge_subsystem* sub
) {
    struct reb_simulation* main_sim = bridge->main_sim;
    struct reb_simulation* sub_sim = sub->sub_sim;
    if (!valid_particle_index(main_sim, sub->host_index)) {
        return bridge_error("host_index is out of bounds during kick");
    }
    const double time_diff = sub_sim->t - main_sim->t;
    const double time_tol = bridge_time_sync_tolerance(main_sim->t, sub_sim->t);
    if (fabs(time_diff) > time_tol) {
        fprintf(
            stderr,
            "[rebound_bridge] main and subsystem times are out of sync: main.t=%.17g sub.t=%.17g diff=%.17g tol=%.17g host_index=%d\n",
            main_sim->t,
            sub_sim->t,
            time_diff,
            time_tol,
            sub->host_index
        );
        return -1;
    }
    /* Snap only the time label inside floating-point tolerance; particle state is unchanged. */
    sub_sim->t = main_sim->t;
    return 0;
}

static int subsystem_total_mass(const struct reb_simulation* sim, double* out) {
    double mt = 0.0;
    for (int i = 0; i < sim->N; i++) {
        mt += sim->particles[i].m;
    }
    if (mt <= 0.0) return bridge_error("subsystem total mass must be positive");
    *out = mt;
    return 0;
}

static int is_other_subsystem_host(
    const struct reb_bridge* bridge,
    const struct reb_bridge_subsystem* current,
    int main_index
) {
    if (bridge->subsystem_coupling != REB_BRIDGE_SUBSYSTEM_COUPLING_RESOLVED) return 0;
    for (int i = 0; i < bridge->n_subsystems; i++) {
        const struct reb_bridge_subsystem* other = &bridge->subsystems[i];
        if (other != current && other->host_index == main_index) return 1;
    }
    return 0;
}

/*
 * Computes differential cross accelerations for one resolved subsystem.
 * By default, other resolved subsystems enter through their host particles'
 * monopoles in main_sim.  In resolved subsystem-coupling mode, those hosts are
 * skipped here and handled by the subsystem-pair residual kick below.
 */
int reb_bridge_internal_compute_subsystem_cross_acc(
    struct reb_bridge* bridge,
    struct reb_bridge_subsystem* sub,
    struct reb_bridge_vec3* internal_acc,
    struct reb_bridge_vec3* main_acc
) {
    struct reb_simulation* main_sim = bridge->main_sim;
    struct reb_simulation* sub_sim = sub->sub_sim;
    const int n = sub_sim->N;
    struct reb_bridge_vec3 bc;
    if (subsystem_barycenter(sub_sim, &bc) != 0) return -1;

    struct reb_bridge_vec3* source_acc = bridge->source_acc;
    memset(internal_acc, 0, sizeof(*internal_acc) * (size_t)n);
    memset(main_acc, 0, sizeof(*main_acc) * (size_t)main_sim->N);

    const struct reb_particle* host = &main_sim->particles[sub->host_index];
    struct reb_bridge_vec3 host_pos = particle_pos(host);
    struct reb_bridge_vec3 bc_global = vec_add(host_pos, bc);
    const double G = main_sim->G;
    double mt = 0.0;
    for (int i = 0; i < n; i++) {
        mt += sub_sim->particles[i].m;
    }
    if (mt <= 0.0) {
        return bridge_error("subsystem total mass must be positive");
    }

    struct reb_bridge_vec3 host_acc = vec3(0.0, 0.0, 0.0);

    /*
     * The bridge kick applies the external tidal field to resolved subsystem
     * members and feeds the equal reaction back to the main-system particles.
     */
    for (int source_index = 0; source_index < main_sim->N; source_index++) {
        if (source_index == sub->host_index) continue;
        if (is_other_subsystem_host(bridge, sub, source_index)) continue;
        const struct reb_particle* source = &main_sim->particles[source_index];
        if (source->m == 0.0) continue;

        struct reb_bridge_vec3 a_mono_on_com = gravity_from_source(bc_global, source, G);
        struct reb_bridge_vec3 a_com_full = vec3(0.0, 0.0, 0.0);
        struct reb_bridge_vec3 source_pos = particle_pos(source);
        struct reb_bridge_vec3 a_source_full = vec3(0.0, 0.0, 0.0);
        for (int i = 0; i < n; i++) {
            const struct reb_particle* member = &sub_sim->particles[i];
            struct reb_bridge_vec3 global_pos =
                vec_add(host_pos, particle_pos(member));
            source_acc[i] = gravity_from_source(global_pos, source, G);
            a_com_full = vec_add(a_com_full, vec_scale(source_acc[i], member->m));
            a_source_full = vec_add(
                a_source_full,
                gravity_from_point_mass(source_pos, global_pos, member->m, G)
            );
        }

        a_com_full = vec_scale(a_com_full, 1.0 / mt);
        host_acc = vec_add(host_acc, vec_sub(a_com_full, a_mono_on_com));

        for (int i = 0; i < n; i++) {
            internal_acc[i] = vec_add(internal_acc[i], vec_sub(source_acc[i], a_com_full));
        }

        struct reb_bridge_vec3 a_source_mono =
            gravity_from_point_mass(source_pos, bc_global, mt, G);
        main_acc[source_index] = vec_add(main_acc[source_index], vec_sub(a_source_full, a_source_mono));
    }

    main_acc[sub->host_index] = vec_add(main_acc[sub->host_index], host_acc);
    return 0;
}

static int apply_resolved_subsystem_pair_cross_kick(
    struct reb_bridge* bridge,
    struct reb_bridge_subsystem* a,
    struct reb_bridge_subsystem* b,
    double dt_kick
) {
    struct reb_simulation* main_sim = bridge->main_sim;
    struct reb_simulation* sub_a = a->sub_sim;
    struct reb_simulation* sub_b = b->sub_sim;

    if (reb_bridge_internal_check_subsystem_time_sync(bridge, a) != 0) return -1;
    if (reb_bridge_internal_check_subsystem_time_sync(bridge, b) != 0) return -1;
    if (!valid_particle_index(main_sim, a->host_index) ||
        !valid_particle_index(main_sim, b->host_index)) {
        return bridge_error("subsystem host_index is out of bounds during resolved coupling");
    }
    if (a->host_index == b->host_index) {
        return bridge_error("resolved subsystem coupling requires distinct host particles");
    }

    const int n_a = sub_a->N;
    const int n_b = sub_b->N;
    const int max_n = n_a > n_b ? n_a : n_b;
    if (reb_bridge_internal_ensure_kick_workspace(bridge, max_n, main_sim->N) != 0) return -1;

    struct reb_bridge_vec3* acc_a = bridge->internal_acc;
    struct reb_bridge_vec3* acc_b = bridge->source_acc;
    memset(acc_a, 0, sizeof(*acc_a) * (size_t)n_a);
    memset(acc_b, 0, sizeof(*acc_b) * (size_t)n_b);

    double mt_a = 0.0;
    double mt_b = 0.0;
    if (subsystem_total_mass(sub_a, &mt_a) != 0) return -1;
    if (subsystem_total_mass(sub_b, &mt_b) != 0) return -1;

    struct reb_bridge_vec3 bc_a;
    struct reb_bridge_vec3 bc_b;
    if (subsystem_barycenter(sub_a, &bc_a) != 0) return -1;
    if (subsystem_barycenter(sub_b, &bc_b) != 0) return -1;

    const struct reb_particle* host_a = &main_sim->particles[a->host_index];
    const struct reb_particle* host_b = &main_sim->particles[b->host_index];
    const struct reb_bridge_vec3 host_pos_a = particle_pos(host_a);
    const struct reb_bridge_vec3 host_pos_b = particle_pos(host_b);
    const struct reb_bridge_vec3 com_a = vec_add(host_pos_a, bc_a);
    const struct reb_bridge_vec3 com_b = vec_add(host_pos_b, bc_b);
    const double G = main_sim->G;

    for (int i = 0; i < n_a; i++) {
        const struct reb_particle* pa = &sub_a->particles[i];
        const struct reb_bridge_vec3 pos_a = vec_add(host_pos_a, particle_pos(pa));
        for (int j = 0; j < n_b; j++) {
            const struct reb_particle* pb = &sub_b->particles[j];
            const struct reb_bridge_vec3 pos_b = vec_add(host_pos_b, particle_pos(pb));
            acc_a[i] = vec_add(
                acc_a[i],
                gravity_from_point_mass(pos_a, pos_b, pb->m, G)
            );
            acc_b[j] = vec_add(
                acc_b[j],
                gravity_from_point_mass(pos_b, pos_a, pa->m, G)
            );
        }
    }

    struct reb_bridge_vec3 com_acc_a = vec3(0.0, 0.0, 0.0);
    struct reb_bridge_vec3 com_acc_b = vec3(0.0, 0.0, 0.0);
    for (int i = 0; i < n_a; i++) {
        com_acc_a = vec_add(com_acc_a, vec_scale(acc_a[i], sub_a->particles[i].m));
    }
    for (int j = 0; j < n_b; j++) {
        com_acc_b = vec_add(com_acc_b, vec_scale(acc_b[j], sub_b->particles[j].m));
    }
    com_acc_a = vec_scale(com_acc_a, 1.0 / mt_a);
    com_acc_b = vec_scale(com_acc_b, 1.0 / mt_b);

    for (int i = 0; i < n_a; i++) {
        struct reb_particle* pa = &sub_a->particles[i];
        const struct reb_bridge_vec3 internal = vec_sub(acc_a[i], com_acc_a);
        pa->vx += internal.x * dt_kick;
        pa->vy += internal.y * dt_kick;
        pa->vz += internal.z * dt_kick;
    }
    for (int j = 0; j < n_b; j++) {
        struct reb_particle* pb = &sub_b->particles[j];
        const struct reb_bridge_vec3 internal = vec_sub(acc_b[j], com_acc_b);
        pb->vx += internal.x * dt_kick;
        pb->vy += internal.y * dt_kick;
        pb->vz += internal.z * dt_kick;
    }

    const struct reb_bridge_vec3 mono_a = gravity_from_source(com_a, host_b, G);
    const struct reb_bridge_vec3 mono_b = gravity_from_source(com_b, host_a, G);
    const struct reb_bridge_vec3 host_acc_a = vec_sub(com_acc_a, mono_a);
    const struct reb_bridge_vec3 host_acc_b = vec_sub(com_acc_b, mono_b);
    struct reb_particle* main_a = &main_sim->particles[a->host_index];
    struct reb_particle* main_b = &main_sim->particles[b->host_index];
    main_a->vx += host_acc_a.x * dt_kick;
    main_a->vy += host_acc_a.y * dt_kick;
    main_a->vz += host_acc_a.z * dt_kick;
    main_b->vx += host_acc_b.x * dt_kick;
    main_b->vy += host_acc_b.y * dt_kick;
    main_b->vz += host_acc_b.z * dt_kick;

    return 0;
}

static int apply_resolved_subsystem_pair_cross_kicks(
    struct reb_bridge* bridge,
    double dt_kick
) {
    if (bridge->subsystem_coupling != REB_BRIDGE_SUBSYSTEM_COUPLING_RESOLVED) return 0;
    for (int i = 0; i < bridge->n_subsystems; i++) {
        for (int j = i + 1; j < bridge->n_subsystems; j++) {
            if (apply_resolved_subsystem_pair_cross_kick(
                    bridge,
                    &bridge->subsystems[i],
                    &bridge->subsystems[j],
                    dt_kick
                ) != 0) {
                return -1;
            }
        }
    }
    for (int i = 0; i < bridge->n_subsystems; i++) {
        reb_simulation_synchronize(bridge->subsystems[i].sub_sim);
    }
    reb_simulation_synchronize(bridge->main_sim);
    return 0;
}

static void apply_subsystem_cross_acc(
    struct reb_bridge* bridge,
    struct reb_bridge_subsystem* sub,
    const struct reb_bridge_vec3* internal_acc,
    const struct reb_bridge_vec3* main_acc,
    double dt_kick
) {
    struct reb_simulation* main_sim = bridge->main_sim;
    struct reb_simulation* sub_sim = sub->sub_sim;
    const int n = sub_sim->N;

    for (int i = 0; i < n; i++) {
        struct reb_particle* p = &sub_sim->particles[i];
        p->vx += internal_acc[i].x * dt_kick;
        p->vy += internal_acc[i].y * dt_kick;
        p->vz += internal_acc[i].z * dt_kick;
    }

    for (int main_index = 0; main_index < main_sim->N; main_index++) {
        struct reb_particle* p = &main_sim->particles[main_index];
        p->vx += main_acc[main_index].x * dt_kick;
        p->vy += main_acc[main_index].y * dt_kick;
        p->vz += main_acc[main_index].z * dt_kick;
    }
}

static int apply_subsystem_cross_kick(
    struct reb_bridge* bridge,
    struct reb_bridge_subsystem* sub,
    double dt_kick,
    double step_dt
) {
    (void)step_dt;
    if (reb_bridge_internal_check_subsystem_time_sync(bridge, sub) != 0) return -1;
    if (reb_bridge_internal_ensure_kick_workspace(bridge, sub->sub_sim->N, bridge->main_sim->N) != 0) return -1;
    if (reb_bridge_internal_compute_subsystem_cross_acc(bridge, sub, bridge->internal_acc, bridge->main_acc) != 0) return -1;
    apply_subsystem_cross_acc(bridge, sub, bridge->internal_acc, bridge->main_acc, dt_kick);

    reb_simulation_synchronize(sub->sub_sim);
    reb_simulation_synchronize(bridge->main_sim);
    return 0;
}

int reb_bridge_apply_cross_kick(struct reb_bridge* bridge, double dt_kick) {
    if (!bridge) return bridge_error("bridge is NULL");
    const double step_dt = 2.0 * dt_kick;
    for (int i = 0; i < bridge->n_subsystems; i++) {
        if (apply_subsystem_cross_kick(bridge, &bridge->subsystems[i], dt_kick, step_dt) != 0) {
            return -1;
        }
    }
    return apply_resolved_subsystem_pair_cross_kicks(bridge, dt_kick);
}

int reb_bridge_internal_apply_cross_kick_for_step(
    struct reb_bridge* bridge,
    double dt_kick,
    double step_dt
) {
    if (!bridge) return bridge_error("bridge is NULL");
    for (int i = 0; i < bridge->n_subsystems; i++) {
        if (apply_subsystem_cross_kick(bridge, &bridge->subsystems[i], dt_kick, step_dt) != 0) {
            return -1;
        }
    }
    return apply_resolved_subsystem_pair_cross_kicks(bridge, dt_kick);
}

