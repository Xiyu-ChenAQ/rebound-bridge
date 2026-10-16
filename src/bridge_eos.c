/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "bridge_internal.h"

#include <float.h>

/* EOS Bridge advance methods. */

static int bridge_advance_by_step_grid(
    struct reb_bridge* bridge,
    double epoch,
    unsigned long long completed,
    unsigned long long steps,
    int (*step)(struct reb_bridge*, double)
) {
    if (!bridge) return bridge_error("bridge is NULL");
    for (unsigned long long i = 0; i < steps; i++) {
        const double target = bridge_outer_grid_time(epoch, completed + i + 1, bridge->dt_outer);
        const double dt = target - bridge->main_sim->t;
        if (!(dt > 0.0)) {
            return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                     "outer step cannot advance time at this epoch; use a relative epoch");
        }
        if (step(bridge, dt) != 0) return -1;
    }
    return 0;
}

int reb_bridge_internal_advance_lf4_2_grid(struct reb_bridge* bridge, double epoch,
                                         unsigned long long completed, unsigned long long count) {
    return bridge_advance_by_step_grid(bridge, epoch, completed, count,
                                      reb_bridge_internal_step_eos_lf4_2);
}

int reb_bridge_internal_step_eos_lf4_2(struct reb_bridge* bridge, double dt) {
    /* a = (1 - 1/sqrt(3))/2 */
    const double a = 0.211324865405187117745425609749;
    const double start = bridge->main_sim->t;
    const double first = start + a * dt;
    const double second = start + (1.0 - a) * dt;
    const double end = start + dt;
    if (dt != 0.0 && !((dt > 0.0 && start < first && first < second && second < end) ||
                      (dt < 0.0 && start > first && first > second && second > end))) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                 "LF4_2 stages cannot be resolved at this epoch; use a relative epoch");
    }

    if (reb_bridge_internal_integrate_all_to_target_symmetric(bridge, first) != 0) return -1;
    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    if (reb_bridge_internal_integrate_all_to_target_symmetric(bridge, second) != 0) return -1;
    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    if (reb_bridge_internal_integrate_all_to_target_symmetric(bridge, end) != 0) return -1;
    return 0;
}

int reb_bridge_internal_step_eos_lf8_6_4(struct reb_bridge* bridge, double dt) {
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
        if (reb_bridge_internal_integrate_all_to_target_symmetric(bridge, bridge->main_sim->t + a[i] * dt) != 0) return -1;
        if (reb_bridge_internal_apply_cross_kick_for_step(bridge, b[i] * dt, dt) != 0) return -1;
    }
    for (int i = 3; i >= 0; i--) {
        if (reb_bridge_internal_integrate_all_to_target_symmetric(bridge, bridge->main_sim->t + a[i] * dt) != 0) return -1;
        if (i > 0) {
            if (reb_bridge_internal_apply_cross_kick_for_step(bridge, b[i - 1] * dt, dt) != 0) return -1;
        }
    }
    return 0;
}

int reb_bridge_advance_eos_lf4_2(struct reb_bridge* bridge, double duration) {
    reb_bridge_clear_error();
    if (!bridge) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "bridge is NULL");
    if (!isfinite(duration) || duration < 0.0 || !isfinite(bridge->main_sim->t + duration)) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                 "duration must be finite and non-negative, with a finite resulting time");
    }
    if (!bridge_duration_on_outer_grid(duration, bridge->dt_outer)) {
        return bridge_error_format(REB_BRIDGE_ERROR_DURATION_NOT_DIVISIBLE,
                                   "duration %.17g must be an integer multiple of dt_outer %.17g",
                                   duration, bridge->dt_outer);
    }
    return reb_bridge_internal_advance_lf4_2_grid(bridge, bridge->main_sim->t, 0,
                                                (unsigned long long)round(duration / bridge->dt_outer));
}

int reb_bridge_advance_eos_lf8_6_4(struct reb_bridge* bridge, double duration) {
    if (!bridge || !bridge_duration_on_outer_grid(duration, bridge->dt_outer)) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "invalid advance duration");
    }
    return bridge_advance_by_step_grid(bridge, bridge->main_sim->t, 0,
                                      (unsigned long long)round(duration / bridge->dt_outer),
                                      reb_bridge_internal_step_eos_lf8_6_4);
}
