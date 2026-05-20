/* SPDX-License-Identifier: GPL-3.0-or-later */

#include "rebound_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct reb_bridge_vec3 {
    double x;
    double y;
    double z;
};

struct reb_bridge_subsystem {
    struct reb_simulation* sub_sim;
    int host_index;
    double dt_inner;
    int owns_sub_sim;
};

struct reb_bridge {
    struct reb_simulation* main_sim;
    double dt_outer;
    int owns_main_sim;
    struct reb_bridge_subsystem* subsystems;
    int n_subsystems;
    int cap_subsystems;
    struct reb_bridge_vec3* internal_acc;
    struct reb_bridge_vec3* source_acc;
    struct reb_bridge_vec3* main_acc;
    int cap_sub_particles;
    int cap_main_particles;
};

static struct reb_bridge_vec3 vec3(double x, double y, double z) {
    struct reb_bridge_vec3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

static struct reb_bridge_vec3 particle_pos(const struct reb_particle* p) {
    return vec3(p->x, p->y, p->z);
}

static struct reb_bridge_vec3 vec_add(struct reb_bridge_vec3 a, struct reb_bridge_vec3 b) {
    return vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static struct reb_bridge_vec3 vec_sub(struct reb_bridge_vec3 a, struct reb_bridge_vec3 b) {
    return vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static struct reb_bridge_vec3 vec_scale(struct reb_bridge_vec3 v, double s) {
    return vec3(v.x * s, v.y * s, v.z * s);
}

static double vec_norm(struct reb_bridge_vec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

static double bridge_time_sync_tolerance(double a, double b) {
    const double scale = fmax(fabs(a), fabs(b));
    return 1.0e-12 + 1.0e-15 * fmax(1.0, scale);
}

static int bridge_error(const char* msg) {
    fprintf(stderr, "[rebound_bridge] %s\n", msg);
    return -1;
}

int reb_bridge_set_integrator_whfast(struct reb_simulation* sim) {
    if (!sim) return bridge_error("sim is NULL");
    struct reb_integrator_whfast_state* whfast =
        (struct reb_integrator_whfast_state*)reb_simulation_set_integrator(sim, "whfast");
    if (!whfast) return bridge_error("failed to set REBOUND integrator to whfast");
    whfast->safe_mode = 1;
    return 0;
}

int reb_bridge_set_integrator_ias15(struct reb_simulation* sim) {
    if (!sim) return bridge_error("sim is NULL");
    struct reb_integrator_ias15_state* ias15 =
        (struct reb_integrator_ias15_state*)reb_simulation_set_integrator(sim, "ias15");
    if (!ias15) return bridge_error("failed to set REBOUND integrator to ias15");
    return 0;
}

static int valid_particle_index(const struct reb_simulation* sim, int index) {
    return sim && index >= 0 && index < sim->N;
}

static struct reb_bridge_body_state body_state(const struct reb_particle* p) {
    struct reb_bridge_body_state s;
    memset(&s, 0, sizeof(s));
    if (!p) return s;
    s.x = p->x;
    s.y = p->y;
    s.z = p->z;
    s.vx = p->vx;
    s.vy = p->vy;
    s.vz = p->vz;
    s.m = p->m;
    return s;
}

static int subsystem_barycenter(
    const struct reb_simulation* sim,
    struct reb_bridge_vec3* out
) {
    if (!sim || !out) return bridge_error("subsystem_barycenter received NULL");
    if (sim->N <= 0) return bridge_error("subsystem has no particles");

    double mt = 0.0;
    struct reb_bridge_vec3 weighted = vec3(0.0, 0.0, 0.0);
    for (int i = 0; i < sim->N; i++) {
        const struct reb_particle* p = &sim->particles[i];
        mt += p->m;
        weighted = vec_add(weighted, vec_scale(particle_pos(p), p->m));
    }
    if (mt <= 0.0) return bridge_error("subsystem total mass must be positive");

    *out = vec_scale(weighted, 1.0 / mt);
    return 0;
}

static struct reb_bridge_vec3 gravity_from_source(
    struct reb_bridge_vec3 pos,
    const struct reb_particle* source,
    double G
) {
    struct reb_bridge_vec3 r = vec_sub(pos, particle_pos(source));
    double rnorm = vec_norm(r);
    double inv_r3 = 1.0 / (rnorm * rnorm * rnorm);
    return vec_scale(r, -G * source->m * inv_r3);
}

static struct reb_bridge_vec3 gravity_from_point_mass(
    struct reb_bridge_vec3 pos,
    struct reb_bridge_vec3 source_pos,
    double source_mass,
    double G
) {
    struct reb_bridge_vec3 r = vec_sub(pos, source_pos);
    double rnorm = vec_norm(r);
    double inv_r3 = 1.0 / (rnorm * rnorm * rnorm);
    return vec_scale(r, -G * source_mass * inv_r3);
}

struct reb_bridge* reb_bridge_create(struct reb_simulation* main_sim, double dt_outer) {
    if (!main_sim) {
        bridge_error("main_sim is NULL");
        return NULL;
    }
    if (dt_outer <= 0.0) {
        bridge_error("dt_outer must be positive");
        return NULL;
    }

    struct reb_bridge* bridge = (struct reb_bridge*)calloc(1, sizeof(*bridge));
    if (!bridge) {
        bridge_error("out of memory");
        return NULL;
    }

    bridge->main_sim = main_sim;
    bridge->dt_outer = dt_outer;
    bridge->owns_main_sim = 0;
    bridge->main_sim->dt = dt_outer;
    return bridge;
}

void reb_bridge_set_owns_main(struct reb_bridge* bridge, int owns_main_sim) {
    if (!bridge) return;
    bridge->owns_main_sim = owns_main_sim ? 1 : 0;
}

void reb_bridge_free(struct reb_bridge* bridge) {
    if (!bridge) return;
    for (int i = 0; i < bridge->n_subsystems; i++) {
        if (bridge->subsystems[i].owns_sub_sim && bridge->subsystems[i].sub_sim) {
            reb_simulation_free(bridge->subsystems[i].sub_sim);
        }
    }
    free(bridge->subsystems);
    free(bridge->internal_acc);
    free(bridge->source_acc);
    free(bridge->main_acc);
    if (bridge->owns_main_sim && bridge->main_sim) {
        reb_simulation_free(bridge->main_sim);
    }
    free(bridge);
}

int reb_bridge_add_subsystem(
    struct reb_bridge* bridge,
    int host_index,
    struct reb_simulation* sub_sim,
    double dt_inner,
    int owns_sub_sim
) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (!sub_sim) return bridge_error("sub_sim is NULL");
    if (!valid_particle_index(bridge->main_sim, host_index)) {
        return bridge_error("host_index is out of bounds");
    }
    if (sub_sim->N <= 0) return bridge_error("sub_sim has no particles");
    if (dt_inner <= 0.0) return bridge_error("dt_inner must be positive");
    if (dt_inner > bridge->dt_outer) return bridge_error("dt_inner must be <= dt_outer");

    double ratio = bridge->dt_outer / dt_inner;
    if (fabs(ratio - round(ratio)) > 1.0e-12) {
        return bridge_error("dt_outer/dt_inner must be an integer ratio");
    }

    if (bridge->n_subsystems == bridge->cap_subsystems) {
        int new_cap = bridge->cap_subsystems == 0 ? 4 : bridge->cap_subsystems * 2;
        struct reb_bridge_subsystem* new_subsystems =
            (struct reb_bridge_subsystem*)realloc(
                bridge->subsystems,
                sizeof(*new_subsystems) * (size_t)new_cap
            );
        if (!new_subsystems) return bridge_error("out of memory");
        bridge->subsystems = new_subsystems;
        bridge->cap_subsystems = new_cap;
    }

    sub_sim->dt = dt_inner;
    struct reb_bridge_subsystem* sub = &bridge->subsystems[bridge->n_subsystems++];
    sub->sub_sim = sub_sim;
    sub->host_index = host_index;
    sub->dt_inner = dt_inner;
    sub->owns_sub_sim = owns_sub_sim ? 1 : 0;
    return 0;
}

int reb_bridge_subsystem_count(const struct reb_bridge* bridge) {
    return bridge ? bridge->n_subsystems : 0;
}

struct reb_simulation* reb_bridge_main_sim(struct reb_bridge* bridge) {
    return bridge ? bridge->main_sim : NULL;
}

struct reb_simulation* reb_bridge_sub_sim(struct reb_bridge* bridge, int subsystem_index) {
    if (!bridge || subsystem_index < 0 || subsystem_index >= bridge->n_subsystems) return NULL;
    return bridge->subsystems[subsystem_index].sub_sim;
}

static int ensure_kick_workspace(struct reb_bridge* bridge, int sub_particles, int main_particles) {
    if (!bridge || sub_particles < 0 || main_particles < 0) {
        return bridge_error("invalid kick workspace request");
    }

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
        bridge->cap_main_particles = main_particles;
    }

    return 0;
}

static int apply_subsystem_cross_kick(
    struct reb_bridge* bridge,
    struct reb_bridge_subsystem* sub,
    double dt_kick
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
    sub_sim->t = main_sim->t;

    const int n = sub_sim->N;
    struct reb_bridge_vec3 bc;
    if (subsystem_barycenter(sub_sim, &bc) != 0) return -1;

    if (ensure_kick_workspace(bridge, n, main_sim->N) != 0) return -1;
    struct reb_bridge_vec3* internal_acc = bridge->internal_acc;
    struct reb_bridge_vec3* source_acc = bridge->source_acc;
    struct reb_bridge_vec3* main_acc = bridge->main_acc;
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

    for (int source_index = 0; source_index < main_sim->N; source_index++) {
        if (source_index == sub->host_index) continue;
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

    for (int i = 0; i < n; i++) {
        struct reb_particle* p = &sub_sim->particles[i];
        p->vx += internal_acc[i].x * dt_kick;
        p->vy += internal_acc[i].y * dt_kick;
        p->vz += internal_acc[i].z * dt_kick;
    }

    struct reb_particle* host_mut = &main_sim->particles[sub->host_index];
    host_mut->vx += host_acc.x * dt_kick;
    host_mut->vy += host_acc.y * dt_kick;
    host_mut->vz += host_acc.z * dt_kick;

    for (int source_index = 0; source_index < main_sim->N; source_index++) {
        if (source_index == sub->host_index) continue;
        struct reb_particle* source = &main_sim->particles[source_index];
        source->vx += main_acc[source_index].x * dt_kick;
        source->vy += main_acc[source_index].y * dt_kick;
        source->vz += main_acc[source_index].z * dt_kick;
    }

    reb_simulation_synchronize(sub_sim);
    reb_simulation_synchronize(main_sim);
    return 0;
}

int reb_bridge_apply_cross_kick(struct reb_bridge* bridge, double dt_kick) {
    if (!bridge) return bridge_error("bridge is NULL");
    for (int i = 0; i < bridge->n_subsystems; i++) {
        if (apply_subsystem_cross_kick(bridge, &bridge->subsystems[i], dt_kick) != 0) {
            return -1;
        }
    }
    return 0;
}

int reb_bridge_step(struct reb_bridge* bridge, double dt) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (dt <= 0.0) return bridge_error("dt must be positive");

    if (reb_bridge_apply_cross_kick(bridge, 0.5 * dt) != 0) return -1;

    const double target = bridge->main_sim->t + dt;
    for (int i = 0; i < bridge->n_subsystems; i++) {
        enum REB_STATUS status = reb_simulation_integrate(bridge->subsystems[i].sub_sim, target);
        if (status != REB_STATUS_SUCCESS) return bridge_error("subsystem integration failed");
    }

    enum REB_STATUS main_status = reb_simulation_integrate(bridge->main_sim, target);
    if (main_status != REB_STATUS_SUCCESS) return bridge_error("main integration failed");

    if (reb_bridge_apply_cross_kick(bridge, 0.5 * dt) != 0) return -1;
    return 0;
}

static int bridge_integrate_all_to_target(struct reb_bridge* bridge, double target) {
    for (int i = 0; i < bridge->n_subsystems; i++) {
        enum REB_STATUS status = reb_simulation_integrate(bridge->subsystems[i].sub_sim, target);
        if (status != REB_STATUS_SUCCESS) return bridge_error("subsystem integration failed");
    }

    enum REB_STATUS main_status = reb_simulation_integrate(bridge->main_sim, target);
    if (main_status != REB_STATUS_SUCCESS) return bridge_error("main integration failed");
    return 0;
}

int reb_bridge_advance(struct reb_bridge* bridge, double duration) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (duration < 0.0) return bridge_error("duration must be non-negative");
    if (duration == 0.0) return 0;

    const double target = bridge->main_sim->t + duration;
    double dt = target - bridge->main_sim->t;
    if (dt > bridge->dt_outer) dt = bridge->dt_outer;

    if (reb_bridge_apply_cross_kick(bridge, 0.5 * dt) != 0) return -1;

    while (1) {
        const double step_target = bridge->main_sim->t + dt;
        if (bridge_integrate_all_to_target(bridge, step_target) != 0) return -1;

        if (bridge->main_sim->t >= target - 1.0e-14) break;

        double next_dt = target - bridge->main_sim->t;
        if (next_dt > bridge->dt_outer) next_dt = bridge->dt_outer;

        if (reb_bridge_apply_cross_kick(bridge, 0.5 * (dt + next_dt)) != 0) return -1;
        dt = next_dt;
    }

    if (reb_bridge_apply_cross_kick(bridge, 0.5 * dt) != 0) return -1;
    return 0;
}

int reb_bridge_snapshot(
    const struct reb_bridge* bridge,
    int subsystem_index,
    struct reb_bridge_snapshot* out
) {
    if (!bridge || !out) return bridge_error("snapshot received NULL");
    if (subsystem_index < 0 || subsystem_index >= bridge->n_subsystems) {
        return bridge_error("subsystem_index is out of bounds");
    }

    const struct reb_bridge_subsystem* sub = &bridge->subsystems[subsystem_index];
    const struct reb_simulation* main_sim = bridge->main_sim;
    const struct reb_simulation* sub_sim = sub->sub_sim;
    if (!valid_particle_index(main_sim, sub->host_index)) {
        return bridge_error("host_index is out of bounds during snapshot");
    }

    memset(out, 0, sizeof(*out));
    out->time = main_sim->t;
    out->host_index = sub->host_index;
    out->member_count = sub_sim->N;
    out->host = body_state(&main_sim->particles[sub->host_index]);

    if (main_sim->N > 0 && sub->host_index != 0) {
        out->host_primary_distance = vec_norm(vec_sub(
            particle_pos(&main_sim->particles[sub->host_index]),
            particle_pos(&main_sim->particles[0])
        ));
    }

    struct reb_bridge_vec3 bc;
    if (subsystem_barycenter(sub_sim, &bc) != 0) return -1;
    out->subsystem_barycenter_distance = vec_norm(bc);

    if (sub_sim->N > 0) {
        out->member0 = body_state(&sub_sim->particles[0]);
    }
    if (sub_sim->N > 1) {
        out->member1 = body_state(&sub_sim->particles[1]);
        out->subsystem_internal_distance = vec_norm(vec_sub(
            particle_pos(&sub_sim->particles[0]),
            particle_pos(&sub_sim->particles[1])
        ));
    }

    return 0;
}

struct reb_simulation* reb_bridge_make_main_sun_emb(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    reb_simulation_add_fmt(sim, "m", 1.0);
    reb_simulation_add_fmt(sim, "m a e", 3.0e-6, 1.0, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

struct reb_simulation* reb_bridge_make_sub_earth_moon(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    const double m_earth = 3.0e-6 * 0.987;
    const double m_moon = 3.0e-6 * 0.013;
    reb_simulation_add_fmt(sim, "m", m_earth);
    reb_simulation_add_fmt(sim, "m a e", m_moon, 0.00257, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

struct reb_bridge* reb_bridge_create_earth_moon(double dt_outer, int sub_ratio) {
    if (sub_ratio <= 0) {
        bridge_error("sub_ratio must be positive");
        return NULL;
    }

    struct reb_simulation* main_sim = reb_bridge_make_main_sun_emb();
    struct reb_simulation* sub_sim = reb_bridge_make_sub_earth_moon();
    if (!main_sim || !sub_sim) {
        if (main_sim) reb_simulation_free(main_sim);
        if (sub_sim) reb_simulation_free(sub_sim);
        return NULL;
    }

    struct reb_bridge* bridge = reb_bridge_create(main_sim, dt_outer);
    if (!bridge) {
        reb_simulation_free(main_sim);
        reb_simulation_free(sub_sim);
        return NULL;
    }

    reb_bridge_set_owns_main(bridge, 1);
    if (reb_bridge_add_subsystem(bridge, 1, sub_sim, dt_outer / (double)sub_ratio, 1) != 0) {
        reb_bridge_free(bridge);
        return NULL;
    }
    return bridge;
}
