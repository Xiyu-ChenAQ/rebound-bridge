/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "bridge_internal.h"

#include <float.h>
#include <limits.h>

/* Order-2 Bridge advance method and method dispatcher. */
static int bridge_duration_is_integral_outer_steps(const struct reb_bridge* bridge, double duration) {
    return bridge_duration_on_outer_grid(duration, bridge->dt_outer);
}
static int bridge_integrator_uses_fixed_steps(const struct reb_simulation* sim) {
    switch (sim->integrator) {
        case REB_INTEGRATOR_WHFAST:
        case REB_INTEGRATOR_SEI:
        case REB_INTEGRATOR_LEAPFROG:
        case REB_INTEGRATOR_NONE:
        case REB_INTEGRATOR_JANUS:
        case REB_INTEGRATOR_MERCURIUS:
        case REB_INTEGRATOR_SABA:
        case REB_INTEGRATOR_EOS:
        case REB_INTEGRATOR_WHFAST512:
        case REB_INTEGRATOR_TRACE:
            return 1;
        case REB_INTEGRATOR_IAS15:
        case REB_INTEGRATOR_BS:
#ifndef REB_BRIDGE_REBOUND_NO_CUSTOM
        case REB_INTEGRATOR_CUSTOM:
#endif
            return 0;
    }
    return 0;
}

static enum REB_STATUS bridge_integrate_fixed_segment(
    struct reb_bridge* bridge,
    struct reb_simulation* sim,
    double target
) {
    const double duration = target - sim->t;
    if (duration == 0.0) return REB_STATUS_SUCCESS;

    const double nominal_dt = fabs(sim->dt);
    if (!(nominal_dt > 0.0) || !isfinite(nominal_dt)) {
        return REB_STATUS_GENERIC_ERROR;
    }

    const double ratio = fabs(duration) / nominal_dt;
    if (!isfinite(ratio)) return REB_STATUS_GENERIC_ERROR;
    const double clock_scale = fmax(fabs(sim->t), fabs(target));
    const double clock_ulp = nextafter(clock_scale, INFINITY) - clock_scale;
    /* Subtracting absolute times loses precision relative to a short segment.
     * Include that rounding uncertainty before ceil: otherwise an integral
     * 40-step segment can become 41 steps at a later epoch, and its reversed
     * counterpart can use a different subdivision. Deliberately non-integral
     * durations outside this rounding bound still round up. */
    const double ratio_tolerance = 64.0 * DBL_EPSILON * fmax(1.0, ratio)
        + 2.0 * clock_ulp / nominal_dt;
    double steps_required = ceil(ratio - ratio_tolerance);
    if (steps_required < 1.0) steps_required = 1.0;
    if (steps_required > (double)UINT_MAX) {
        return REB_STATUS_GENERIC_ERROR;
    }
    const unsigned int steps = (unsigned int)steps_required;

    const double saved_dt = sim->dt;
    const int saved_exact_finish_time = sim->exact_finish_time;
    /* 阶段内均分步长，不在末尾补短步 */
    const double segment_dt = duration / (double)steps;
    const double stop_time = target - 0.5 * segment_dt;
    /* REBOUND stops on its floating-point clock. Bound cumulative rounding
     * of its two half-clock updates so that it cannot cross the midpoint
     * stopping boundary one step early or late. Near that limit, rejecting
     * the epoch is safer than relabelling an incorrectly evolved state. */
    if (sim->t + 0.5 * segment_dt == sim->t || stop_time == target ||
        2.0 * (double)steps * clock_ulp >= 0.25 * fabs(segment_dt)) {
        bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                          "component substeps cannot be counted reliably at this epoch; use a relative epoch");
        return REB_STATUS_GENERIC_ERROR;
    }
    sim->dt = segment_dt;
    sim->exact_finish_time = 0;
    const unsigned long long previous_steps = sim->steps_done;
    const enum REB_STATUS status = reb_bridge_internal_integrate_component(bridge, sim, stop_time);
    sim->dt = saved_dt;
    sim->exact_finish_time = saved_exact_finish_time;

    if (status != REB_STATUS_SUCCESS) {
        fprintf(
            stderr,
            "[rebound_bridge] fixed segment failed: status=%d t=%.17g target=%.17g "
            "duration=%.17g steps=%u segment_dt=%.17g\n",
            (int)status,
            sim->t,
            target,
            duration,
            steps,
            segment_dt
        );
        return status;
    }
    if (sim->steps_done - previous_steps != (unsigned long long)steps) {
        bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                          "component executed an unexpected number of substeps; use a relative epoch");
        return REB_STATUS_GENERIC_ERROR;
    }
    const double time_scale = fmax(1.0, fmax(fabs(sim->t), fabs(target)));
    const double time_tolerance = fmax(
        bridge_time_advance_tolerance(sim->t, target),
        256.0 * DBL_EPSILON * (double)steps * time_scale
    );
    if (fabs(sim->t - target) > time_tolerance) {
        fprintf(
            stderr,
            "[rebound_bridge] fixed segment missed target: t=%.17g target=%.17g "
            "diff=%.17g duration=%.17g steps=%u segment_dt=%.17g\n",
            sim->t,
            target,
            sim->t - target,
            duration,
            steps,
            segment_dt
        );
        return REB_STATUS_GENERIC_ERROR;
    }
    sim->t = target;
    return REB_STATUS_SUCCESS;
}

static enum REB_STATUS bridge_integrate_to_target(
    struct reb_bridge* bridge,
    struct reb_simulation* sim,
    double target
) {
    if (bridge_integrator_uses_fixed_steps(sim)) {
        return bridge_integrate_fixed_segment(bridge, sim, target);
    }
#ifndef REB_BRIDGE_PYTHON
    if (target != sim->t &&
        sim->t + copysign(fabs(sim->dt), target - sim->t) == sim->t) {
        bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                          "component timestep cannot resolve time at this epoch; use a relative epoch");
        return REB_STATUS_GENERIC_ERROR;
    }
#endif
    return reb_bridge_internal_integrate_component(bridge, sim, target);
}

static int bridge_component_integration_error(
    const struct reb_bridge* bridge, int subsystem_index, enum REB_STATUS status
) {
    /* Preserve specific Bridge failures, including clock-resolution checks. */
    if (reb_bridge_last_error() != REB_BRIDGE_ERROR_NONE) return -1;

    const char* reason;
    switch (status) {
        case REB_STATUS_GENERIC_ERROR: reason = "generic integration error"; break;
        case REB_STATUS_NO_PARTICLES: reason = "no particles remain"; break;
        case REB_STATUS_ENCOUNTER: reason = "close encounter (exit_min_distance)"; break;
        case REB_STATUS_ESCAPE: reason = "particle escaped (exit_max_distance)"; break;
        case REB_STATUS_USER: reason = "user requested stop"; break;
        case REB_STATUS_SIGINT: reason = "integration interrupted (SIGINT)"; break;
        case REB_STATUS_COLLISION: reason = "particle collision"; break;
        default: reason = "unexpected integration status"; break;
    }
    if (subsystem_index < 0) {
        return bridge_error_format(REB_BRIDGE_ERROR_INTERNAL,
                                   "main integration failed: REBOUND status %d (%s)",
                                   (int)status, reason);
    }
    return bridge_error_format(REB_BRIDGE_ERROR_INTERNAL,
                               "subsystem %d (host particle %d) integration failed: REBOUND status %d (%s)",
                               subsystem_index, bridge->subsystems[subsystem_index].host_index,
                               (int)status, reason);
}

int reb_bridge_internal_integrate_all_to_target(struct reb_bridge* bridge, double target) {
    for (int i = 0; i < bridge->n_subsystems; i++) {
        enum REB_STATUS status = reb_bridge_internal_integrate_component(bridge, bridge->subsystems[i].sub_sim, target);
        if (status != REB_STATUS_SUCCESS) return bridge_component_integration_error(bridge, i, status);
    }

    enum REB_STATUS main_status = reb_bridge_internal_integrate_component(bridge, bridge->main_sim, target);
    if (main_status != REB_STATUS_SUCCESS) return bridge_component_integration_error(bridge, -1, main_status);
    return 0;
}

int reb_bridge_internal_integrate_all_to_target_symmetric(
    struct reb_bridge* bridge,
    double target
) {
    for (int i = 0; i < bridge->n_subsystems; i++) {
        enum REB_STATUS status = bridge_integrate_to_target(bridge, bridge->subsystems[i].sub_sim, target);
        if (status != REB_STATUS_SUCCESS) {
            return bridge_component_integration_error(bridge, i, status);
        }
    }

    enum REB_STATUS main_status = bridge_integrate_to_target(bridge, bridge->main_sim, target);
    if (main_status != REB_STATUS_SUCCESS) {
        return bridge_component_integration_error(bridge, -1, main_status);
    }
    return 0;
}

void reb_bridge_internal_snap_all_times(struct reb_bridge* bridge, double time) {
    bridge->main_sim->t = time;
    for (int i = 0; i < bridge->n_subsystems; i++) {
        bridge->subsystems[i].sub_sim->t = time;
    }
}

int reb_bridge_internal_step_order2(struct reb_bridge* bridge, double dt) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    if (reb_bridge_internal_integrate_all_to_target_symmetric(bridge, bridge->main_sim->t + dt) != 0) return -1;
    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    return 0;
}

int reb_bridge_step(struct reb_bridge* bridge, double dt) {
    reb_bridge_clear_error();
    if (!bridge) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "bridge is NULL");
    if (!isfinite(dt) || !isfinite(bridge->main_sim->t + dt)) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "dt and resulting time must be finite");
    }
    if (dt == 0.0) return 0;
    if (bridge->main_sim->t + dt == bridge->main_sim->t) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                 "step cannot advance time at this epoch; use a relative epoch");
    }
    switch (bridge->method) {
        case REB_BRIDGE_METHOD_ORDER2:
            return reb_bridge_internal_step_order2(bridge, dt);
        case REB_BRIDGE_METHOD_EOS_LF4_2:
            return reb_bridge_internal_step_eos_lf4_2(bridge, dt);
        case REB_BRIDGE_METHOD_INVALID:
            return bridge_error("invalid bridge method");
    }
    return bridge_error("invalid bridge method");
}

static int bridge_advance_order2_grid(struct reb_bridge* bridge, double epoch,
                                      unsigned long long completed, unsigned long long steps) {
    if (steps == 0) return 0;
    double dt = bridge_outer_grid_time(epoch, completed + 1, bridge->dt_outer) - bridge->main_sim->t;
    if (!(dt > 0.0)) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                 "outer step cannot advance time at this epoch; use a relative epoch");
    }

    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;

    /* Adjacent plain half-kicks are merged, preserving the KDK composition over many outer steps. */
    for (unsigned long long i = 0; i < steps; i++) {
        const double step_target = bridge_outer_grid_time(epoch, completed + i + 1, bridge->dt_outer);
        if (reb_bridge_internal_integrate_all_to_target_symmetric(bridge, step_target) != 0) return -1;

        if (i + 1 == steps) break;

        const double next_dt = bridge_outer_grid_time(epoch, completed + i + 2, bridge->dt_outer) - step_target;
        if (!(next_dt > 0.0)) {
            return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                     "outer grid is not distinguishable at this epoch; use a relative epoch");
        }

        if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * (dt + next_dt), 0.5 * (dt + next_dt) * 2.0) != 0) return -1;
        dt = next_dt;
    }

    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    return 0;
}

int reb_bridge_internal_advance_grid(struct reb_bridge* bridge, double epoch,
                                   unsigned long long completed, unsigned long long count) {
    reb_bridge_clear_error();
    /* 步数要能被 double 精确表示 */
    if (!bridge || !isfinite(epoch) || completed > 9007199254740991ULL ||
        count > 9007199254740991ULL - completed ||
        !isfinite(bridge_outer_grid_time(epoch, completed + count, bridge->dt_outer))) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "invalid absolute outer-step grid");
    }
    switch (bridge->method) {
        case REB_BRIDGE_METHOD_ORDER2:
            return bridge_advance_order2_grid(bridge, epoch, completed, count);
        case REB_BRIDGE_METHOD_EOS_LF4_2:
            return reb_bridge_internal_advance_lf4_2_grid(bridge, epoch, completed, count);
        case REB_BRIDGE_METHOD_INVALID:
            return bridge_error("invalid bridge method");
    }
    return bridge_error("invalid bridge method");
}

int reb_bridge_advance(struct reb_bridge* bridge, double duration) {
    reb_bridge_clear_error();
    if (!bridge) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "bridge is NULL");
    if (!isfinite(duration) || duration < 0.0 || !isfinite(bridge->main_sim->t + duration)) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                 "duration must be finite and non-negative, with a finite resulting time");
    }
    if (!bridge_duration_is_integral_outer_steps(bridge, duration)) {
        return bridge_error_format(REB_BRIDGE_ERROR_DURATION_NOT_DIVISIBLE,
                                   "duration %.17g must be an integer multiple of dt_outer %.17g",
                                   duration, bridge->dt_outer);
    }
    return reb_bridge_internal_advance_grid(bridge, bridge->main_sim->t, 0,
                                           (unsigned long long)round(duration / bridge->dt_outer));
}
