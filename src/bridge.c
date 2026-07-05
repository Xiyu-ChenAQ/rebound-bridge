/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "bridge_internal.h"

/* Core Bridge lifecycle and stable API. */

int reb_bridge_set_integrator_whfast(struct reb_simulation* sim) {
    if (!sim) return bridge_error("sim is NULL");
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->ri_whfast.safe_mode = 1;
    return 0;
}

int reb_bridge_set_integrator_ias15(struct reb_simulation* sim) {
    if (!sim) return bridge_error("sim is NULL");
    sim->integrator = REB_INTEGRATOR_IAS15;
    return 0;
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
    bridge->method = REB_BRIDGE_METHOD_ORDER2;
    bridge->cross_kick_kernel = REB_BRIDGE_CROSS_KICK_PLAIN;
    bridge->force_gradient_scale = 1.0;
    bridge->corrector = REB_BRIDGE_CORRECTOR_NONE;
    bridge->corrector_scale = 1.0;
    bridge->main_sim->dt = dt_outer;
    return bridge;
}

void reb_bridge_set_owns_main(struct reb_bridge* bridge, int owns_main_sim) {
    if (!bridge) return;
    bridge->owns_main_sim = owns_main_sim ? 1 : 0;
}

static int valid_bridge_method(enum reb_bridge_method method) {
    switch (method) {
        case REB_BRIDGE_METHOD_ORDER2:
        case REB_BRIDGE_METHOD_EOS_LF4_2:
        case REB_BRIDGE_METHOD_EOS_LF8_6_4:
        case REB_BRIDGE_METHOD_FORWARD_FG4:
            return 1;
        case REB_BRIDGE_METHOD_INVALID:
            return 0;
    }
    return 0;
}

int reb_bridge_set_method(struct reb_bridge* bridge, enum reb_bridge_method method) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (!valid_bridge_method(method)) return bridge_error("invalid bridge method");
    bridge->method = method;
    return 0;
}

enum reb_bridge_method reb_bridge_get_method(const struct reb_bridge* bridge) {
    return bridge ? bridge->method : REB_BRIDGE_METHOD_INVALID;
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
    free(bridge->fg_internal_acc);
    free(bridge->fg_main_acc);
    free(bridge->fd_internal_acc);
    free(bridge->fd_main_acc);
    free(bridge->saved_sub_pos);
    free(bridge->saved_main_pos);
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
