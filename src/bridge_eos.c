/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "bridge_internal.h"

/* EOS Bridge advance methods. */

static int bridge_advance_by_step(
    struct reb_bridge* bridge,
    double duration,
    int (*step)(struct reb_bridge*, double)
) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (duration < 0.0) return bridge_error("duration must be non-negative");
    if (duration == 0.0) return 0;

    const double target = bridge->main_sim->t + duration;
    while (!bridge_reached_time(bridge->main_sim->t, target)) {
        double dt = target - bridge->main_sim->t;
        if (dt > bridge->dt_outer) dt = bridge->dt_outer;
        if (step(bridge, dt) != 0) return -1;
    }
    return 0;
}

static int bridge_eos_lf4_2_step(struct reb_bridge* bridge, double dt) {
    const double a = 0.211324865405187117745425609749;

    if (reb_bridge_internal_integrate_all_to_target(bridge, bridge->main_sim->t + a * dt) != 0) return -1;
    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    if (reb_bridge_internal_integrate_all_to_target(bridge, bridge->main_sim->t + (1.0 - 2.0 * a) * dt) != 0) return -1;
    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    if (reb_bridge_internal_integrate_all_to_target(bridge, bridge->main_sim->t + a * dt) != 0) return -1;
    return 0;
}

static int bridge_eos_lf8_6_4_step(struct reb_bridge* bridge, double dt) {
    static const double a[4] = {
        0.0711334264982231177779387300061549964174,
        0.241153427956640098736487795326289649618,
        0.521411761772814789212136078067994229991,
        -0.333698616227678005726562603400438876027
    };
    static const double b[4] = {
        0.183083687472197221961703757166430291072,
        0.310782859898574869507522291054262796375,
        -0.0265646185119588006972121379164987592663,
        0.0653961422823734184559721793911134363710
    };

    for (int i = 0; i < 4; i++) {
        if (reb_bridge_internal_integrate_all_to_target(bridge, bridge->main_sim->t + a[i] * dt) != 0) return -1;
        if (reb_bridge_internal_apply_cross_kick_for_step(bridge, b[i] * dt, dt) != 0) return -1;
    }
    for (int i = 3; i >= 0; i--) {
        if (reb_bridge_internal_integrate_all_to_target(bridge, bridge->main_sim->t + a[i] * dt) != 0) return -1;
        if (i > 0) {
            if (reb_bridge_internal_apply_cross_kick_for_step(bridge, b[i - 1] * dt, dt) != 0) return -1;
        }
    }
    return 0;
}

int reb_bridge_advance_eos_lf4_2(struct reb_bridge* bridge, double duration) {
    return bridge_advance_by_step(bridge, duration, bridge_eos_lf4_2_step);
}

int reb_bridge_advance_eos_lf8_6_4(struct reb_bridge* bridge, double duration) {
    return bridge_advance_by_step(bridge, duration, bridge_eos_lf8_6_4_step);
}
