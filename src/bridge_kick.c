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

/*
 * Computes differential cross accelerations for one resolved subsystem.
 * Other resolved subsystems enter only through their host particles' monopoles
 * in main_sim; there is no direct resolved-resolved kick between subsystems.
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
    double step_dt,
    enum reb_bridge_cross_kick_kernel kernel
) {
    if (reb_bridge_internal_check_subsystem_time_sync(bridge, sub) != 0) return -1;
    if (reb_bridge_internal_ensure_kick_workspace(bridge, sub->sub_sim->N, bridge->main_sim->N) != 0) return -1;
    if (reb_bridge_internal_compute_subsystem_cross_acc(bridge, sub, bridge->internal_acc, bridge->main_acc) != 0) return -1;
    if (kernel == REB_BRIDGE_CROSS_KICK_FORCE_GRADIENT) {
        if (reb_bridge_internal_apply_force_gradient_correction(
                bridge,
                sub,
                bridge->internal_acc,
                bridge->main_acc,
                step_dt
            ) != 0) {
            return -1;
        }
    }
    apply_subsystem_cross_acc(bridge, sub, bridge->internal_acc, bridge->main_acc, dt_kick);

    reb_simulation_synchronize(sub->sub_sim);
    reb_simulation_synchronize(bridge->main_sim);
    return 0;
}

int reb_bridge_apply_cross_kick(struct reb_bridge* bridge, double dt_kick) {
    if (!bridge) return bridge_error("bridge is NULL");
    const double step_dt = 2.0 * dt_kick;
    for (int i = 0; i < bridge->n_subsystems; i++) {
        if (apply_subsystem_cross_kick(bridge, &bridge->subsystems[i], dt_kick, step_dt, bridge->cross_kick_kernel) != 0) {
            return -1;
        }
    }
    return 0;
}

int reb_bridge_internal_apply_cross_kick_for_step(
    struct reb_bridge* bridge,
    double dt_kick,
    double step_dt
) {
    if (!bridge) return bridge_error("bridge is NULL");
    for (int i = 0; i < bridge->n_subsystems; i++) {
        if (apply_subsystem_cross_kick(bridge, &bridge->subsystems[i], dt_kick, step_dt, bridge->cross_kick_kernel) != 0) {
            return -1;
        }
    }
    return 0;
}

int reb_bridge_internal_apply_cross_kick_with_kernel_for_step(
    struct reb_bridge* bridge,
    double dt_kick,
    double step_dt,
    enum reb_bridge_cross_kick_kernel kernel
) {
    if (!bridge) return bridge_error("bridge is NULL");
    for (int i = 0; i < bridge->n_subsystems; i++) {
        if (apply_subsystem_cross_kick(bridge, &bridge->subsystems[i], dt_kick, step_dt, kernel) != 0) {
            return -1;
        }
    }
    return 0;
}
