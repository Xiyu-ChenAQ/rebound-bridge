/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "bridge_internal.h"

/* Experimental Bridge-level commutator processor. */

static int valid_corrector(enum reb_bridge_corrector corrector) {
    switch (corrector) {
        case REB_BRIDGE_CORRECTOR_NONE:
        case REB_BRIDGE_CORRECTOR_COMMUTATOR_AB:
            return 1;
        case REB_BRIDGE_CORRECTOR_INVALID:
            return 0;
    }
    return 0;
}

int reb_bridge_set_corrector(struct reb_bridge* bridge, enum reb_bridge_corrector corrector) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (!valid_corrector(corrector)) return bridge_error("invalid Bridge-level corrector");
    bridge->corrector = corrector;
    return 0;
}

enum reb_bridge_corrector reb_bridge_get_corrector(const struct reb_bridge* bridge) {
    return bridge ? bridge->corrector : REB_BRIDGE_CORRECTOR_INVALID;
}

int reb_bridge_set_corrector_scale(struct reb_bridge* bridge, double scale) {
    if (!bridge) return bridge_error("bridge is NULL");
    if (!isfinite(scale)) return bridge_error("Bridge-level corrector scale must be finite");
    bridge->corrector_scale = scale;
    return 0;
}

double reb_bridge_get_corrector_scale(const struct reb_bridge* bridge) {
    return bridge ? bridge->corrector_scale : NAN;
}

static int bridge_corrector_A(struct reb_bridge* bridge, double dt) {
    if (dt == 0.0) return 0;
    return reb_bridge_internal_integrate_all_to_target(bridge, bridge->main_sim->t + dt);
}

static int bridge_corrector_B(struct reb_bridge* bridge, double dt_kick) {
    if (dt_kick == 0.0) return 0;
    return reb_bridge_internal_apply_cross_kick_for_step(bridge, dt_kick, 2.0 * fabs(dt_kick));
}

static int bridge_apply_commutator_ab_corrector(
    struct reb_bridge* bridge,
    enum reb_bridge_corrector_direction direction
) {
    if (bridge->corrector_scale == 0.0) return 0;

    const double t0 = bridge->main_sim->t;
    const double magnitude = sqrt(fabs(bridge->corrector_scale)) * bridge->dt_outer;
    const double a = magnitude;
    const double b = copysign(magnitude, bridge->corrector_scale);
    int status = 0;

    if (direction == REB_BRIDGE_CORRECTOR_FORWARD) {
        status = bridge_corrector_A(bridge, a);
        if (status == 0) status = bridge_corrector_B(bridge, b);
        if (status == 0) status = bridge_corrector_A(bridge, -a);
        if (status == 0) status = bridge_corrector_B(bridge, -b);
    } else if (direction == REB_BRIDGE_CORRECTOR_INVERSE) {
        status = bridge_corrector_B(bridge, b);
        if (status == 0) status = bridge_corrector_A(bridge, a);
        if (status == 0) status = bridge_corrector_B(bridge, -b);
        if (status == 0) status = bridge_corrector_A(bridge, -a);
    } else {
        return bridge_error("invalid Bridge-level corrector direction");
    }

    if (status != 0) return status;
    reb_bridge_internal_snap_all_times(bridge, t0);
    return 0;
}

int reb_bridge_apply_corrector(
    struct reb_bridge* bridge,
    enum reb_bridge_corrector_direction direction
) {
    if (!bridge) return bridge_error("bridge is NULL");
    switch (bridge->corrector) {
        case REB_BRIDGE_CORRECTOR_NONE:
            return 0;
        case REB_BRIDGE_CORRECTOR_COMMUTATOR_AB:
            return bridge_apply_commutator_ab_corrector(bridge, direction);
        case REB_BRIDGE_CORRECTOR_INVALID:
            return bridge_error("invalid Bridge-level corrector");
    }
    return bridge_error("invalid Bridge-level corrector");
}
