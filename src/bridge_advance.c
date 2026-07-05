/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "bridge_internal.h"

/* Order-2 Bridge advance method and method dispatcher. */
int reb_bridge_internal_integrate_all_to_target(struct reb_bridge* bridge, double target) {
    for (int i = 0; i < bridge->n_subsystems; i++) {
        enum REB_STATUS status = reb_simulation_integrate(bridge->subsystems[i].sub_sim, target);
        if (status != REB_STATUS_SUCCESS) return bridge_error("subsystem integration failed");
    }

    enum REB_STATUS main_status = reb_simulation_integrate(bridge->main_sim, target);
    if (main_status != REB_STATUS_SUCCESS) return bridge_error("main integration failed");
    return 0;
}

void reb_bridge_internal_snap_all_times(struct reb_bridge* bridge, double time) {
    bridge->main_sim->t = time;
    for (int i = 0; i < bridge->n_subsystems; i++) {
        bridge->subsystems[i].sub_sim->t = time;
    }
}

int reb_bridge_step(struct reb_bridge* bridge, double dt) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    if (reb_bridge_internal_integrate_all_to_target(bridge, bridge->main_sim->t + dt) != 0) return -1;
    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    return 0;
}

static int bridge_advance_order2(struct reb_bridge* bridge, double duration) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (duration < 0.0) return bridge_error("duration must be non-negative");
    if (duration == 0.0) return 0;

    const double target = bridge->main_sim->t + duration;
    double dt = target - bridge->main_sim->t;
    if (dt > bridge->dt_outer) dt = bridge->dt_outer;

    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;

    /* Adjacent plain half-kicks are merged, preserving the KDK composition over many outer steps. */
    while (1) {
        const double step_target = bridge->main_sim->t + dt;
        if (reb_bridge_internal_integrate_all_to_target(bridge, step_target) != 0) return -1;

        if (bridge_reached_time(bridge->main_sim->t, target)) break;

        double next_dt = target - bridge->main_sim->t;
        if (next_dt > bridge->dt_outer) next_dt = bridge->dt_outer;

        if (bridge->cross_kick_kernel == REB_BRIDGE_CROSS_KICK_FORCE_GRADIENT) {
            if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
            if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * next_dt, next_dt) != 0) return -1;
        } else {
            if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * (dt + next_dt), 0.5 * (dt + next_dt) * 2.0) != 0) return -1;
        }
        dt = next_dt;
    }

    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    return 0;
}

int reb_bridge_advance(struct reb_bridge* bridge, double duration) {
    if (!bridge) return bridge_error("bridge is NULL");
    switch (bridge->method) {
        case REB_BRIDGE_METHOD_ORDER2:
            return bridge_advance_order2(bridge, duration);
        case REB_BRIDGE_METHOD_EOS_LF4_2:
            return reb_bridge_advance_eos_lf4_2(bridge, duration);
        case REB_BRIDGE_METHOD_EOS_LF8_6_4:
            return reb_bridge_advance_eos_lf8_6_4(bridge, duration);
        case REB_BRIDGE_METHOD_FORWARD_FG4:
            return reb_bridge_advance_forward_fg4(bridge, duration);
        case REB_BRIDGE_METHOD_INVALID:
            return bridge_error("invalid bridge method");
    }
    return bridge_error("invalid bridge method");
}
