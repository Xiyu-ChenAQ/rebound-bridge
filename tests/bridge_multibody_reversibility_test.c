/* SPDX-License-Identifier: GPL-3.0-or-later */

#include "rebound_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct test_state {
    struct reb_bridge* bridge;
    struct reb_simulation* main_sim;
    struct reb_simulation* sub_sim;
};

static struct test_state make_state(enum reb_bridge_method method) {
    struct test_state state = {0};
    const double dt_outer = 4.0 / 365.25;
    const double m_jupiter = 9.547919384243266e-4;
    const double m_io = 4.704e-8;
    const double m_europa = 2.528e-8;
    const double m_ganymede = 7.804e-8;
    const double subsystem_mass = m_jupiter + m_io + m_europa + m_ganymede;

    state.main_sim = reb_simulation_create();
    state.sub_sim = reb_simulation_create();
    if (!state.main_sim || !state.sub_sim) return state;
    state.main_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    state.sub_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(state.main_sim) != 0 ||
        reb_bridge_set_integrator_whfast(state.sub_sim) != 0) {
        return state;
    }

    reb_simulation_add_fmt(state.main_sim, "m", subsystem_mass);
    reb_simulation_add_fmt(state.sub_sim, "m", m_jupiter);
    reb_simulation_add_fmt(state.sub_sim, "m a e f", m_io, 0.002819, 0.004, 0.2);
    reb_simulation_add_fmt(state.sub_sim, "m a e f", m_europa, 0.004486, 0.009, 2.0);
    reb_simulation_add_fmt(state.sub_sim, "m a e f", m_ganymede, 0.007155, 0.001, 4.0);
    reb_simulation_move_to_com(state.sub_sim);

    state.bridge = reb_bridge_create(state.main_sim, dt_outer);
    if (!state.bridge) return state;
    reb_bridge_set_owns_main(state.bridge, 1);
    if (reb_bridge_set_method(state.bridge, method) != 0 ||
        reb_bridge_add_subsystem(state.bridge, 0, state.sub_sim, dt_outer / 40.0, 1) != 0) {
        return state;
    }
    return state;
}

static void free_state(struct test_state* state) {
    if (state->bridge) {
        reb_bridge_free(state->bridge);
    } else {
        if (state->main_sim) reb_simulation_free(state->main_sim);
        if (state->sub_sim) reb_simulation_free(state->sub_sim);
    }
}

static double simulation_delta(
    const struct reb_simulation* simulation,
    const struct reb_simulation* reference
) {
    double maximum = fabs(simulation->t - reference->t);
    if (simulation->N != reference->N) return INFINITY;
    for (int i = 0; i < simulation->N; i++) {
        const struct reb_particle* a = &simulation->particles[i];
        const struct reb_particle* b = &reference->particles[i];
        const double av[] = {a->x, a->y, a->z, a->vx, a->vy, a->vz};
        const double bv[] = {b->x, b->y, b->z, b->vx, b->vy, b->vz};
        for (int j = 0; j < 6; j++) {
            const double delta = fabs(av[j] - bv[j]);
            if (delta > maximum) maximum = delta;
        }
    }
    return maximum;
}

static int measure_roundtrip(enum reb_bridge_method method, double* delta_out) {
    const double dt_outer = 4.0 / 365.25;
    const int cycles = 10;
    struct test_state state = make_state(method);
    if (!state.bridge) {
        free_state(&state);
        return -1;
    }

    reb_simulation_synchronize(state.main_sim);
    reb_simulation_synchronize(state.sub_sim);
    struct reb_simulation* initial_main = reb_simulation_copy(state.main_sim);
    struct reb_simulation* initial_sub = reb_simulation_copy(state.sub_sim);
    if (!initial_main || !initial_sub) {
        if (initial_main) reb_simulation_free(initial_main);
        if (initial_sub) reb_simulation_free(initial_sub);
        free_state(&state);
        return -1;
    }

    for (int i = 0; i < cycles; i++) {
        if (reb_bridge_step(state.bridge, dt_outer) != 0) return -1;
    }
    for (int i = 0; i < cycles; i++) {
        if (reb_bridge_step(state.bridge, -dt_outer) != 0) return -1;
    }
    reb_simulation_synchronize(state.main_sim);
    reb_simulation_synchronize(state.sub_sim);
    *delta_out = fmax(
        simulation_delta(state.main_sim, initial_main),
        simulation_delta(state.sub_sim, initial_sub)
    );

    reb_simulation_free(initial_main);
    reb_simulation_free(initial_sub);
    free_state(&state);
    return 0;
}

int main(void) {
    double order2_delta = 0.0;
    double lf42_delta = 0.0;
    if (measure_roundtrip(REB_BRIDGE_METHOD_ORDER2, &order2_delta) != 0 ||
        measure_roundtrip(REB_BRIDGE_METHOD_EOS_LF4_2, &lf42_delta) != 0) {
        fprintf(stderr, "failed to run multibody reversibility test\n");
        return EXIT_FAILURE;
    }

    printf("multibody_order2_roundtrip_max=%.17g\n", order2_delta);
    printf("multibody_lf42_roundtrip_max=%.17g\n", lf42_delta);
    if (order2_delta > 1.0e-9 || lf42_delta > 1.0e-9) {
        fprintf(
            stderr,
            "multibody roundtrip drift too large: order2=%.17g lf42=%.17g\n",
            order2_delta,
            lf42_delta
        );
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
