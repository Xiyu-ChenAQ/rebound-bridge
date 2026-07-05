/* SPDX-License-Identifier: GPL-3.0-or-later */

#include "rebound_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int fail(const char* msg) {
    fprintf(stderr, "%s\n", msg);
    return EXIT_FAILURE;
}

static double particle_state_delta(const struct reb_particle* a, const struct reb_particle* b) {
    double max_delta = 0.0;
    const double values_a[] = {a->x, a->y, a->z, a->vx, a->vy, a->vz};
    const double values_b[] = {b->x, b->y, b->z, b->vx, b->vy, b->vz};
    for (int i = 0; i < 6; i++) {
        const double delta = fabs(values_a[i] - values_b[i]);
        if (delta > max_delta) max_delta = delta;
    }
    return max_delta;
}

static struct reb_bridge* make_bridge(struct reb_simulation** main_sim_out) {
    struct reb_simulation* sim = reb_simulation_create();

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    struct reb_particle sun = {0};
    sun.m = 1.0;
    reb_simulation_add(sim, sun);

    struct reb_bridge* bridge = reb_bridge_create(sim, 0.01);
    if (!bridge) {
        reb_simulation_free(sim);
        return NULL;
    }
    reb_bridge_set_owns_main(bridge, 1);
    *main_sim_out = sim;
    return bridge;
}

static struct reb_simulation* make_two_body_subsystem(void) {
    struct reb_simulation* sub = reb_simulation_create();
    sub->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sub) != 0) {
        reb_simulation_free(sub);
        return NULL;
    }

    struct reb_particle primary = {0};
    primary.m = 3.0e-6;
    reb_simulation_add(sub, primary);

    struct reb_particle secondary = {0};
    secondary.m = 3.7e-8;
    secondary.x = 0.0025;
    secondary.vy = 0.22;
    reb_simulation_add(sub, secondary);

    reb_simulation_move_to_com(sub);
    return sub;
}

int main(void) {
    if (reb_bridge_get_method(NULL) != REB_BRIDGE_METHOD_INVALID) {
        return fail("NULL bridge should report REB_BRIDGE_METHOD_INVALID");
    }

    struct reb_simulation* main_sim = NULL;
    struct reb_bridge* bridge = make_bridge(&main_sim);
    if (!bridge || !main_sim) {
        return fail("failed to create test bridge");
    }

    if (reb_bridge_get_method(bridge) != REB_BRIDGE_METHOD_ORDER2) {
        reb_bridge_free(bridge);
        return fail("new bridge should default to REB_BRIDGE_METHOD_ORDER2");
    }
    if (reb_bridge_get_cross_kick_kernel(bridge) != REB_BRIDGE_CROSS_KICK_PLAIN) {
        reb_bridge_free(bridge);
        return fail("new bridge should default to plain cross-kick kernel");
    }
    if (reb_bridge_get_corrector(bridge) != REB_BRIDGE_CORRECTOR_NONE) {
        reb_bridge_free(bridge);
        return fail("new bridge should default to no Bridge-level corrector");
    }
    if (fabs(reb_bridge_get_force_gradient_scale(bridge) - 1.0) > 0.0) {
        reb_bridge_free(bridge);
        return fail("new bridge should default to unit force-gradient scale");
    }
    if (fabs(reb_bridge_get_corrector_scale(bridge) - 1.0) > 0.0) {
        reb_bridge_free(bridge);
        return fail("new bridge should default to unit corrector scale");
    }

    if (reb_bridge_set_method(bridge, REB_BRIDGE_METHOD_EOS_LF4_2) != 0) {
        reb_bridge_free(bridge);
        return fail("setting EOS LF4_2 method should succeed");
    }
    if (reb_bridge_get_method(bridge) != REB_BRIDGE_METHOD_EOS_LF4_2) {
        reb_bridge_free(bridge);
        return fail("bridge should remember EOS LF4_2 method");
    }

    if (reb_bridge_set_method(bridge, REB_BRIDGE_METHOD_EOS_LF8_6_4) != 0) {
        reb_bridge_free(bridge);
        return fail("setting EOS LF8_6_4 method should succeed");
    }
    if (reb_bridge_get_method(bridge) != REB_BRIDGE_METHOD_EOS_LF8_6_4) {
        reb_bridge_free(bridge);
        return fail("bridge should remember EOS LF8_6_4 method");
    }

    if (reb_bridge_set_method(bridge, REB_BRIDGE_METHOD_FORWARD_FG4) != 0) {
        reb_bridge_free(bridge);
        return fail("setting forward force-gradient fourth-order method should succeed");
    }
    if (reb_bridge_get_method(bridge) != REB_BRIDGE_METHOD_FORWARD_FG4) {
        reb_bridge_free(bridge);
        return fail("bridge should remember forward force-gradient fourth-order method");
    }

    if (reb_bridge_set_method(bridge, (enum reb_bridge_method)999) == 0) {
        reb_bridge_free(bridge);
        return fail("invalid bridge method should fail");
    }
    if (reb_bridge_get_method(bridge) != REB_BRIDGE_METHOD_FORWARD_FG4) {
        reb_bridge_free(bridge);
        return fail("invalid bridge method should not change current method");
    }
    if (reb_bridge_set_cross_kick_kernel(bridge, REB_BRIDGE_CROSS_KICK_FORCE_GRADIENT) != 0) {
        reb_bridge_free(bridge);
        return fail("setting force-gradient cross-kick kernel should succeed");
    }
    if (reb_bridge_get_cross_kick_kernel(bridge) != REB_BRIDGE_CROSS_KICK_FORCE_GRADIENT) {
        reb_bridge_free(bridge);
        return fail("bridge should remember force-gradient cross-kick kernel");
    }
    if (reb_bridge_set_cross_kick_kernel(bridge, (enum reb_bridge_cross_kick_kernel)999) == 0) {
        reb_bridge_free(bridge);
        return fail("invalid cross-kick kernel should fail");
    }
    if (reb_bridge_get_cross_kick_kernel(bridge) != REB_BRIDGE_CROSS_KICK_FORCE_GRADIENT) {
        reb_bridge_free(bridge);
        return fail("invalid cross-kick kernel should not change current kernel");
    }
    if (reb_bridge_set_corrector(bridge, REB_BRIDGE_CORRECTOR_COMMUTATOR_AB) != 0) {
        reb_bridge_free(bridge);
        return fail("setting commutator Bridge-level corrector should succeed");
    }
    if (reb_bridge_get_corrector(bridge) != REB_BRIDGE_CORRECTOR_COMMUTATOR_AB) {
        reb_bridge_free(bridge);
        return fail("bridge should remember commutator Bridge-level corrector");
    }
    if (reb_bridge_set_corrector(bridge, (enum reb_bridge_corrector)999) == 0) {
        reb_bridge_free(bridge);
        return fail("invalid Bridge-level corrector should fail");
    }
    if (reb_bridge_get_corrector(bridge) != REB_BRIDGE_CORRECTOR_COMMUTATOR_AB) {
        reb_bridge_free(bridge);
        return fail("invalid Bridge-level corrector should not change current corrector");
    }
    if (reb_bridge_set_corrector_scale(bridge, -0.25) != 0) {
        reb_bridge_free(bridge);
        return fail("finite Bridge-level corrector scale should be accepted");
    }
    if (fabs(reb_bridge_get_corrector_scale(bridge) + 0.25) > 1.0e-15) {
        reb_bridge_free(bridge);
        return fail("bridge should remember Bridge-level corrector scale");
    }
    if (reb_bridge_set_force_gradient_scale(bridge, 0.25) != 0) {
        reb_bridge_free(bridge);
        return fail("finite force-gradient scale should be accepted");
    }
    if (fabs(reb_bridge_get_force_gradient_scale(bridge) - 0.25) > 1.0e-15) {
        reb_bridge_free(bridge);
        return fail("bridge should remember force-gradient scale");
    }

    struct reb_particle host = {0};
    host.m = 3.037e-6;
    host.x = 1.0;
    host.vy = 2.0 * M_PI;
    reb_simulation_add(main_sim, host);

    struct reb_simulation* sub = make_two_body_subsystem();
    if (!sub) {
        reb_bridge_free(bridge);
        return fail("failed to create test subsystem");
    }
    if (reb_bridge_add_subsystem(bridge, 1, sub, 0.01, 1) != 0) {
        reb_simulation_free(sub);
        reb_bridge_free(bridge);
        return fail("failed to add test subsystem");
    }
    if (reb_bridge_set_cross_kick_kernel(bridge, REB_BRIDGE_CROSS_KICK_PLAIN) != 0 ||
        reb_bridge_set_corrector_scale(bridge, 0.25) != 0 ||
        reb_bridge_apply_corrector(bridge, REB_BRIDGE_CORRECTOR_FORWARD) != 0) {
        reb_bridge_free(bridge);
        return fail("Bridge-level corrector warm-start should succeed");
    }
    const struct reb_particle main_after_warm = main_sim->particles[1];
    const struct reb_particle sub0_after_warm = sub->particles[0];
    if (reb_bridge_apply_corrector(bridge, REB_BRIDGE_CORRECTOR_INVERSE) != 0 ||
        reb_bridge_apply_corrector(bridge, REB_BRIDGE_CORRECTOR_FORWARD) != 0) {
        reb_bridge_free(bridge);
        return fail("Bridge-level corrector inverse/forward cycle should succeed");
    }
    if (particle_state_delta(&main_after_warm, &main_sim->particles[1]) > 1.0e-10 ||
        particle_state_delta(&sub0_after_warm, &sub->particles[0]) > 1.0e-10) {
        reb_bridge_free(bridge);
        return fail("Bridge-level corrector inverse/forward cycle should return to processed state");
    }
    if (reb_bridge_apply_corrector(bridge, REB_BRIDGE_CORRECTOR_INVERSE) != 0) {
        reb_bridge_free(bridge);
        return fail("Bridge-level corrector cold output should succeed");
    }

    struct reb_bridge_force_gradient_check fg_check;
    if (reb_bridge_force_gradient_finite_difference_check(bridge, 0, 1.0e-7, &fg_check) != 0) {
        reb_bridge_free(bridge);
        return fail("force-gradient finite-difference check should succeed");
    }
    if (!isfinite(fg_check.max_abs_difference) || !isfinite(fg_check.max_relative_difference)) {
        reb_bridge_free(bridge);
        return fail("force-gradient finite-difference check should return finite metrics");
    }
    if (!isfinite(fg_check.internal_weighted_acc_residual) ||
        !isfinite(fg_check.main_weighted_acc_residual) ||
        !isfinite(fg_check.force_gradient_internal_weighted_acc_residual) ||
        !isfinite(fg_check.force_gradient_main_weighted_acc_residual)) {
        reb_bridge_free(bridge);
        return fail("force-gradient check should return finite canonical residual metrics");
    }
    if (fg_check.internal_weighted_acc_residual > 1.0e-20 ||
        fg_check.force_gradient_internal_weighted_acc_residual > 1.0e-15) {
        reb_bridge_free(bridge);
        return fail("subsystem cross-kick and force-gradient should preserve internal COM momentum");
    }

    if (reb_bridge_advance(bridge, 0.02) != 0) {
        reb_bridge_free(bridge);
        return fail("configured reb_bridge_advance should succeed");
    }
    if (fabs(main_sim->t - 0.02) > 1.0e-12) {
        reb_bridge_free(bridge);
        return fail("configured reb_bridge_advance should advance to requested time");
    }

    reb_bridge_free(bridge);
    return EXIT_SUCCESS;
}
