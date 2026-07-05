/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "bridge_internal.h"

/* Forward fourth-order force-gradient Bridge advance method. */

static int bridge_forward_fg4_step(struct reb_bridge* bridge, double dt) {
    if (reb_bridge_internal_apply_cross_kick_with_kernel_for_step(
            bridge,
            dt / 6.0,
            dt,
            REB_BRIDGE_CROSS_KICK_PLAIN
        ) != 0) {
        return -1;
    }
    if (reb_bridge_internal_integrate_all_to_target(bridge, bridge->main_sim->t + 0.5 * dt) != 0) return -1;
    if (reb_bridge_internal_apply_cross_kick_with_kernel_for_step(
            bridge,
            (2.0 / 3.0) * dt,
            dt,
            REB_BRIDGE_CROSS_KICK_FORCE_GRADIENT
        ) != 0) {
        return -1;
    }
    if (reb_bridge_internal_integrate_all_to_target(bridge, bridge->main_sim->t + 0.5 * dt) != 0) return -1;
    if (reb_bridge_internal_apply_cross_kick_with_kernel_for_step(
            bridge,
            dt / 6.0,
            dt,
            REB_BRIDGE_CROSS_KICK_PLAIN
        ) != 0) {
        return -1;
    }
    return 0;
}

int reb_bridge_advance_forward_fg4(struct reb_bridge* bridge, double duration) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (duration < 0.0) return bridge_error("duration must be non-negative");
    if (duration == 0.0) return 0;

    const double target = bridge->main_sim->t + duration;
    while (!bridge_reached_time(bridge->main_sim->t, target)) {
        double dt = target - bridge->main_sim->t;
        if (dt > bridge->dt_outer) dt = bridge->dt_outer;
        if (bridge_forward_fg4_step(bridge, dt) != 0) return -1;
    }
    return 0;
}
