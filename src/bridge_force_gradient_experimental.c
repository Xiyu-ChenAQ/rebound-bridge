/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "bridge_internal.h"

/* Experimental force-gradient cross-kick kernel and diagnostics. */

static double mass_weighted_acc_residual(
    const struct reb_simulation* sim,
    const struct reb_bridge_vec3* acc
) {
    struct reb_bridge_vec3 weighted = vec3(0.0, 0.0, 0.0);
    if (!sim || !acc) return INFINITY;
    for (int i = 0; i < sim->N; i++) {
        weighted = vec_add(weighted, vec_scale(acc[i], sim->particles[i].m));
    }
    return vec_norm(weighted);
}

static int valid_cross_kick_kernel(enum reb_bridge_cross_kick_kernel kernel) {
    switch (kernel) {
        case REB_BRIDGE_CROSS_KICK_PLAIN:
        case REB_BRIDGE_CROSS_KICK_FORCE_GRADIENT:
            return 1;
        case REB_BRIDGE_CROSS_KICK_INVALID:
            return 0;
    }
    return 0;
}

int reb_bridge_set_cross_kick_kernel(
    struct reb_bridge* bridge,
    enum reb_bridge_cross_kick_kernel kernel
) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (!valid_cross_kick_kernel(kernel)) return bridge_error("invalid cross-kick kernel");
    bridge->cross_kick_kernel = kernel;
    return 0;
}

enum reb_bridge_cross_kick_kernel reb_bridge_get_cross_kick_kernel(const struct reb_bridge* bridge) {
    return bridge ? bridge->cross_kick_kernel : REB_BRIDGE_CROSS_KICK_INVALID;
}

int reb_bridge_set_force_gradient_scale(struct reb_bridge* bridge, double scale) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (!isfinite(scale)) return bridge_error("force-gradient scale must be finite");
    bridge->force_gradient_scale = scale;
    return 0;
}

double reb_bridge_get_force_gradient_scale(const struct reb_bridge* bridge) {
    return bridge ? bridge->force_gradient_scale : NAN;
}

static int bridge_force_gradient_coefficient(const struct reb_bridge* bridge, double* coefficient) {
    if (!bridge || !coefficient) return bridge_error("force-gradient coefficient received NULL");
    switch (bridge->method) {
        case REB_BRIDGE_METHOD_ORDER2:
            *coefficient = -1.0 / 12.0;
            return 0;
        case REB_BRIDGE_METHOD_EOS_LF4_2:
            *coefficient = (2.0 - sqrt(3.0)) / 12.0;
            return 0;
        case REB_BRIDGE_METHOD_EOS_LF8_6_4:
            return bridge_error("force-gradient cross-kick is not defined for EOS LF8_6_4");
        case REB_BRIDGE_METHOD_FORWARD_FG4:
            *coefficient = 1.0 / 24.0;
            return 0;
        case REB_BRIDGE_METHOD_INVALID:
            return bridge_error("invalid bridge method");
    }
    return bridge_error("invalid bridge method");
}

static struct reb_bridge_vec3 tidal_tensor_apply(
    struct reb_bridge_vec3 r,
    struct reb_bridge_vec3 dr
) {
    const double r2 = vec_dot(r, r);
    const double rinv = 1.0 / sqrt(r2);
    const double rinv3 = rinv * rinv * rinv;
    const double rinv5 = rinv3 * rinv * rinv;
    const double projection = 3.0 * vec_dot(r, dr) * rinv5;
    return vec_sub(vec_scale(dr, rinv3), vec_scale(r, projection));
}

static int compute_subsystem_cross_acc_directional_derivative(
    struct reb_bridge* bridge,
    struct reb_bridge_subsystem* sub,
    const struct reb_bridge_vec3* internal_acc,
    const struct reb_bridge_vec3* main_acc,
    struct reb_bridge_vec3* derivative_internal_acc,
    struct reb_bridge_vec3* derivative_main_acc
) {
    struct reb_simulation* main_sim = bridge->main_sim;
    struct reb_simulation* sub_sim = sub->sub_sim;
    const int n = sub_sim->N;
    struct reb_bridge_vec3 bc;
    if (subsystem_barycenter(sub_sim, &bc) != 0) return -1;

    memset(derivative_internal_acc, 0, sizeof(*derivative_internal_acc) * (size_t)n);
    memset(derivative_main_acc, 0, sizeof(*derivative_main_acc) * (size_t)main_sim->N);

    const struct reb_particle* host = &main_sim->particles[sub->host_index];
    const struct reb_bridge_vec3 host_pos = particle_pos(host);
    const struct reb_bridge_vec3 bc_global = vec_add(host_pos, bc);
    const struct reb_bridge_vec3 host_b = main_acc[sub->host_index];
    const double G = main_sim->G;
    double mt = 0.0;
    for (int i = 0; i < n; i++) {
        mt += sub_sim->particles[i].m;
    }
    if (mt <= 0.0) {
        return bridge_error("subsystem total mass must be positive");
    }

    for (int source_index = 0; source_index < main_sim->N; source_index++) {
        if (source_index == sub->host_index) continue;
        const struct reb_particle* source = &main_sim->particles[source_index];
        if (source->m == 0.0) continue;

        const struct reb_bridge_vec3 source_pos = particle_pos(source);
        const struct reb_bridge_vec3 source_b = main_acc[source_index];
        const struct reb_bridge_vec3 dR = vec_sub(host_b, source_b);
        const struct reb_bridge_vec3 R = vec_sub(bc_global, source_pos);
        const struct reb_bridge_vec3 tidal_mono = tidal_tensor_apply(R, dR);
        const struct reb_bridge_vec3 d_a_com_mono = vec_scale(tidal_mono, -G * source->m);
        const struct reb_bridge_vec3 d_a_source_mono = vec_scale(tidal_mono, G * mt);
        struct reb_bridge_vec3 d_a_com_full = vec3(0.0, 0.0, 0.0);
        struct reb_bridge_vec3 d_a_source_full = vec3(0.0, 0.0, 0.0);

        for (int i = 0; i < n; i++) {
            const struct reb_particle* member = &sub_sim->particles[i];
            const struct reb_bridge_vec3 member_pos = vec_add(host_pos, particle_pos(member));
            const struct reb_bridge_vec3 r = vec_sub(member_pos, source_pos);
            const struct reb_bridge_vec3 member_b_global = vec_add(host_b, internal_acc[i]);
            const struct reb_bridge_vec3 dr = vec_sub(member_b_global, source_b);
            const struct reb_bridge_vec3 tidal_full = tidal_tensor_apply(r, dr);
            const struct reb_bridge_vec3 d_a_member = vec_scale(tidal_full, -G * source->m);
            bridge->source_acc[i] = d_a_member;
            d_a_com_full = vec_add(d_a_com_full, vec_scale(d_a_member, member->m));
            d_a_source_full = vec_add(
                d_a_source_full,
                vec_scale(tidal_full, G * member->m)
            );
        }

        d_a_com_full = vec_scale(d_a_com_full, 1.0 / mt);
        derivative_main_acc[sub->host_index] = vec_add(
            derivative_main_acc[sub->host_index],
            vec_sub(d_a_com_full, d_a_com_mono)
        );
        derivative_main_acc[source_index] = vec_add(
            derivative_main_acc[source_index],
            vec_sub(d_a_source_full, d_a_source_mono)
        );

        for (int i = 0; i < n; i++) {
            derivative_internal_acc[i] = vec_add(
                derivative_internal_acc[i],
                vec_sub(bridge->source_acc[i], d_a_com_full)
            );
        }
    }
    return 0;
}

int reb_bridge_internal_apply_force_gradient_correction(
    struct reb_bridge* bridge,
    struct reb_bridge_subsystem* sub,
    struct reb_bridge_vec3* internal_acc,
    struct reb_bridge_vec3* main_acc,
    double step_dt
) {
    struct reb_simulation* main_sim = bridge->main_sim;
    struct reb_simulation* sub_sim = sub->sub_sim;
    const int n = sub_sim->N;
    if (compute_subsystem_cross_acc_directional_derivative(
            bridge,
            sub,
            internal_acc,
            main_acc,
            bridge->fg_internal_acc,
            bridge->fg_main_acc
        ) != 0) {
        return -1;
    }

    double coefficient = 0.0;
    if (bridge_force_gradient_coefficient(bridge, &coefficient) != 0) return -1;
    const double correction = coefficient * bridge->force_gradient_scale * step_dt * step_dt;
    for (int i = 0; i < n; i++) {
        internal_acc[i] = vec_add(internal_acc[i], vec_scale(bridge->fg_internal_acc[i], correction));
    }
    for (int i = 0; i < main_sim->N; i++) {
        main_acc[i] = vec_add(main_acc[i], vec_scale(bridge->fg_main_acc[i], correction));
    }
    return 0;
}

static void update_force_gradient_check_metric(
    struct reb_bridge_force_gradient_check* out,
    struct reb_bridge_vec3 analytic,
    struct reb_bridge_vec3 finite_difference
) {
    const struct reb_bridge_vec3 diff = vec_sub(analytic, finite_difference);
    const double abs_diff = vec_norm(diff);
    const double abs_analytic = vec_norm(analytic);
    const double abs_fd = vec_norm(finite_difference);
    const double denom = fmax(1.0e-300, fmax(abs_analytic, abs_fd));
    const double rel_diff = abs_diff / denom;
    if (abs_diff > out->max_abs_difference) out->max_abs_difference = abs_diff;
    if (rel_diff > out->max_relative_difference) out->max_relative_difference = rel_diff;
    if (abs_analytic > out->max_abs_analytic) out->max_abs_analytic = abs_analytic;
    if (abs_fd > out->max_abs_finite_difference) out->max_abs_finite_difference = abs_fd;
}

int reb_bridge_force_gradient_finite_difference_check(
    struct reb_bridge* bridge,
    int subsystem_index,
    double epsilon,
    struct reb_bridge_force_gradient_check* out
) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (!out) return bridge_error("force-gradient check output is NULL");
    if (subsystem_index < 0 || subsystem_index >= bridge->n_subsystems) {
        return bridge_error("subsystem_index is out of bounds");
    }
    if (epsilon <= 0.0) return bridge_error("epsilon must be positive");

    struct reb_bridge_subsystem* sub = &bridge->subsystems[subsystem_index];
    struct reb_simulation* sub_sim = sub->sub_sim;
    struct reb_simulation* main_sim = bridge->main_sim;
    const int n = sub_sim->N;
    if (reb_bridge_internal_check_subsystem_time_sync(bridge, sub) != 0) return -1;
    if (reb_bridge_internal_ensure_kick_workspace(bridge, n, main_sim->N) != 0) return -1;

    if (reb_bridge_internal_compute_subsystem_cross_acc(bridge, sub, bridge->internal_acc, bridge->main_acc) != 0) return -1;
    if (compute_subsystem_cross_acc_directional_derivative(
            bridge,
            sub,
            bridge->internal_acc,
            bridge->main_acc,
            bridge->fg_internal_acc,
            bridge->fg_main_acc
        ) != 0) {
        return -1;
    }

    for (int i = 0; i < n; i++) {
        bridge->saved_sub_pos[i] = particle_pos(&sub_sim->particles[i]);
        sub_sim->particles[i].x += epsilon * bridge->internal_acc[i].x;
        sub_sim->particles[i].y += epsilon * bridge->internal_acc[i].y;
        sub_sim->particles[i].z += epsilon * bridge->internal_acc[i].z;
    }
    for (int i = 0; i < main_sim->N; i++) {
        bridge->saved_main_pos[i] = particle_pos(&main_sim->particles[i]);
        main_sim->particles[i].x += epsilon * bridge->main_acc[i].x;
        main_sim->particles[i].y += epsilon * bridge->main_acc[i].y;
        main_sim->particles[i].z += epsilon * bridge->main_acc[i].z;
    }

    int status = reb_bridge_internal_compute_subsystem_cross_acc(bridge, sub, bridge->fd_internal_acc, bridge->fd_main_acc);

    for (int i = 0; i < n; i++) {
        sub_sim->particles[i].x = bridge->saved_sub_pos[i].x;
        sub_sim->particles[i].y = bridge->saved_sub_pos[i].y;
        sub_sim->particles[i].z = bridge->saved_sub_pos[i].z;
    }
    for (int i = 0; i < main_sim->N; i++) {
        main_sim->particles[i].x = bridge->saved_main_pos[i].x;
        main_sim->particles[i].y = bridge->saved_main_pos[i].y;
        main_sim->particles[i].z = bridge->saved_main_pos[i].z;
    }
    if (status != 0) return -1;

    memset(out, 0, sizeof(*out));
    out->internal_weighted_acc_residual = mass_weighted_acc_residual(sub_sim, bridge->internal_acc);
    out->main_weighted_acc_residual = mass_weighted_acc_residual(main_sim, bridge->main_acc);
    out->force_gradient_internal_weighted_acc_residual = mass_weighted_acc_residual(sub_sim, bridge->fg_internal_acc);
    out->force_gradient_main_weighted_acc_residual = mass_weighted_acc_residual(main_sim, bridge->fg_main_acc);
    for (int i = 0; i < n; i++) {
        const struct reb_bridge_vec3 fd = vec_scale(
            vec_sub(bridge->fd_internal_acc[i], bridge->internal_acc[i]),
            1.0 / epsilon
        );
        update_force_gradient_check_metric(out, bridge->fg_internal_acc[i], fd);
    }
    for (int i = 0; i < main_sim->N; i++) {
        const struct reb_bridge_vec3 fd = vec_scale(
            vec_sub(bridge->fd_main_acc[i], bridge->main_acc[i]),
            1.0 / epsilon
        );
        update_force_gradient_check_metric(out, bridge->fg_main_acc[i], fd);
    }
    return 0;
}
