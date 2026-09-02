/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "bridge_internal.h"

#include <float.h>
#include <limits.h>

/* Order-2 Bridge advance method and method dispatcher. */
static int bridge_duration_is_integral_outer_steps(const struct reb_bridge* bridge, double duration) {
    if (!isfinite(duration) || duration < 0.0) return 0;
    const double steps = duration / bridge->dt_outer;
    return isfinite(steps) && fabs(steps - round(steps)) <= 64.0 * DBL_EPSILON * fmax(1.0, steps);
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
        case REB_INTEGRATOR_CUSTOM:
            return 0;
    }
    return 0;
}

static enum REB_STATUS bridge_integrate_fixed_segment(
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
    const double ratio_tolerance = 64.0 * DBL_EPSILON * fmax(1.0, ratio);
    double steps_required = ceil(ratio - ratio_tolerance);
    if (steps_required < 1.0) steps_required = 1.0;
    if (steps_required > (double)UINT_MAX) {
        return REB_STATUS_GENERIC_ERROR;
    }
    const unsigned int steps = (unsigned int)steps_required;

    const double saved_dt = sim->dt;
    const int saved_exact_finish_time = sim->exact_finish_time;
    const double segment_dt = duration / (double)steps;
    const double stop_time = target - 0.5 * segment_dt;
    sim->dt = segment_dt;
    sim->exact_finish_time = 0;
    const enum REB_STATUS status = reb_simulation_integrate(sim, stop_time);
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
    struct reb_simulation* sim,
    double target
) {
    if (bridge_integrator_uses_fixed_steps(sim)) {
        return bridge_integrate_fixed_segment(sim, target);
    }
    return reb_simulation_integrate(sim, target);
}

int reb_bridge_internal_integrate_all_to_target(struct reb_bridge* bridge, double target) {
    for (int i = 0; i < bridge->n_subsystems; i++) {
        enum REB_STATUS status = reb_simulation_integrate(bridge->subsystems[i].sub_sim, target);
        if (status != REB_STATUS_SUCCESS) return bridge_error("subsystem integration failed");
    }

    enum REB_STATUS main_status = reb_simulation_integrate(bridge->main_sim, target);
    if (main_status != REB_STATUS_SUCCESS) return bridge_error("main integration failed");
    return 0;
}

int reb_bridge_internal_integrate_all_to_target_symmetric(
    struct reb_bridge* bridge,
    double target
) {
    for (int i = 0; i < bridge->n_subsystems; i++) {
        enum REB_STATUS status = bridge_integrate_to_target(bridge->subsystems[i].sub_sim, target);
        if (status != REB_STATUS_SUCCESS) return bridge_error("subsystem integration failed");
    }

    enum REB_STATUS main_status = bridge_integrate_to_target(bridge->main_sim, target);
    if (main_status != REB_STATUS_SUCCESS) return bridge_error("main integration failed");
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
    if (reb_bridge_internal_integrate_all_to_target(bridge, bridge->main_sim->t + dt) != 0) return -1;
    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    return 0;
}

int reb_bridge_step(struct reb_bridge* bridge, double dt) {
    if (!bridge) return bridge_error("bridge is NULL");
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

static int bridge_advance_order2(struct reb_bridge* bridge, double duration) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (!bridge_duration_is_integral_outer_steps(bridge, duration)) {
        bridge_error("duration must be an integer multiple of dt_outer");
        reb_bridge_error_state = REB_BRIDGE_ERROR_DURATION_NOT_DIVISIBLE;
        return -1;
    }
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

        if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * (dt + next_dt), 0.5 * (dt + next_dt) * 2.0) != 0) return -1;
        dt = next_dt;
    }

    if (reb_bridge_internal_apply_cross_kick_for_step(bridge, 0.5 * dt, dt) != 0) return -1;
    return 0;
}

int reb_bridge_advance(struct reb_bridge* bridge, double duration) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (!bridge_duration_is_integral_outer_steps(bridge, duration)) {
        bridge_error("duration must be an integer multiple of dt_outer");
        reb_bridge_error_state = REB_BRIDGE_ERROR_DURATION_NOT_DIVISIBLE;
        return -1;
    }
    switch (bridge->method) {
        case REB_BRIDGE_METHOD_ORDER2:
            return bridge_advance_order2(bridge, duration);
        case REB_BRIDGE_METHOD_EOS_LF4_2:
            return reb_bridge_advance_eos_lf4_2(bridge, duration);
        case REB_BRIDGE_METHOD_INVALID:
            return bridge_error("invalid bridge method");
    }
    return bridge_error("invalid bridge method");
}
