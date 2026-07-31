/* SPDX-License-Identifier: GPL-3.0-or-later */

#include "rebound_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct bridge_state {
    struct reb_bridge* bridge;
    struct reb_simulation* main_sim;
    struct reb_simulation* sub_sim;
};

static int fail(const char* message) {
    fprintf(stderr, "%s\n", message);
    return EXIT_FAILURE;
}

static struct bridge_state make_state(enum reb_bridge_method method) {
    struct bridge_state state = {0};
    const double dt_outer = 0.01;
    const double dt_inner = dt_outer / 8.0;
    const double subsystem_mass = 3.037e-6;

    state.main_sim = reb_simulation_create();
    state.sub_sim = reb_simulation_create();
    if (!state.main_sim || !state.sub_sim) return state;

    state.main_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    state.sub_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(state.main_sim) != 0 ||
        reb_bridge_set_integrator_whfast(state.sub_sim) != 0) {
        return state;
    }

    reb_simulation_add_fmt(state.main_sim, "m", 1.0);
    reb_simulation_add_fmt(state.main_sim, "m a e", subsystem_mass, 1.0, 0.05);
    reb_simulation_move_to_com(state.main_sim);

    reb_simulation_add_fmt(state.sub_sim, "m", 3.0e-6);
    reb_simulation_add_fmt(state.sub_sim, "m a e", 3.7e-8, 0.0025, 0.1);
    reb_simulation_move_to_com(state.sub_sim);

    state.bridge = reb_bridge_create(state.main_sim, dt_outer);
    if (!state.bridge) return state;
    reb_bridge_set_owns_main(state.bridge, 1);
    if (reb_bridge_set_method(state.bridge, method) != 0 ||
        reb_bridge_add_subsystem(state.bridge, 1, state.sub_sim, dt_inner, 1) != 0) {
        return state;
    }
    return state;
}

static void free_state(struct bridge_state* state) {
    if (state->bridge) {
        reb_bridge_free(state->bridge);
    } else {
        if (state->main_sim) reb_simulation_free(state->main_sim);
        if (state->sub_sim) reb_simulation_free(state->sub_sim);
    }
    *state = (struct bridge_state){0};
}

static double simulation_delta(
    const struct reb_simulation* a,
    const struct reb_simulation* b
) {
    if (!a || !b || a->N != b->N) return INFINITY;
    double max_delta = fabs(a->t - b->t);
    for (int i = 0; i < a->N; i++) {
        const struct reb_particle* pa = &a->particles[i];
        const struct reb_particle* pb = &b->particles[i];
        const double av[] = {pa->x, pa->y, pa->z, pa->vx, pa->vy, pa->vz};
        const double bv[] = {pb->x, pb->y, pb->z, pb->vx, pb->vy, pb->vz};
        for (int j = 0; j < 6; j++) {
            const double delta = fabs(av[j] - bv[j]);
            if (delta > max_delta) max_delta = delta;
        }
    }
    return max_delta;
}

int main(void) {
    const double dt_outer = 0.01;
    struct bridge_state order2 = make_state(REB_BRIDGE_METHOD_ORDER2);
    struct bridge_state lf_dispatch = make_state(REB_BRIDGE_METHOD_EOS_LF4_2);
    struct bridge_state lf_direct = make_state(REB_BRIDGE_METHOD_EOS_LF4_2);
    struct bridge_state lf_tiny_step = make_state(REB_BRIDGE_METHOD_EOS_LF4_2);
    if (!order2.bridge || !lf_dispatch.bridge || !lf_direct.bridge || !lf_tiny_step.bridge) {
        free_state(&order2);
        free_state(&lf_dispatch);
        free_state(&lf_direct);
        free_state(&lf_tiny_step);
        return fail("failed to create step-dispatch test states");
    }

    if (reb_bridge_step(order2.bridge, dt_outer) != 0 ||
        reb_bridge_step(lf_dispatch.bridge, dt_outer) != 0 ||
        reb_bridge_advance_eos_lf4_2(lf_direct.bridge, dt_outer) != 0 ||
        reb_bridge_step(lf_tiny_step.bridge, 1.0e-20) != 0) {
        free_state(&order2);
        free_state(&lf_dispatch);
        free_state(&lf_direct);
        free_state(&lf_tiny_step);
        return fail("step-dispatch test integration failed");
    }

    const double method_delta = fmax(
        simulation_delta(order2.main_sim, lf_dispatch.main_sim),
        simulation_delta(order2.sub_sim, lf_dispatch.sub_sim)
    );
    const double dispatch_delta = fmax(
        simulation_delta(lf_dispatch.main_sim, lf_direct.main_sim),
        simulation_delta(lf_dispatch.sub_sim, lf_direct.sub_sim)
    );

    free_state(&order2);
    free_state(&lf_dispatch);
    free_state(&lf_direct);
    free_state(&lf_tiny_step);

    if (method_delta <= 1.0e-15) {
        return fail("reb_bridge_step ignored the configured Bridge method");
    }
    if (dispatch_delta > 1.0e-14) {
        fprintf(stderr, "LF step dispatch differs from direct LF advance: %.17g\n", dispatch_delta);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
