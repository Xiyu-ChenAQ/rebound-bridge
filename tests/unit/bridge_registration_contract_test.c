/* SPDX-License-Identifier: GPL-3.0-or-later */

#include "rebound_bridge.h"

#include <stdio.h>
#include <stdlib.h>

enum contract_case {
    CONTRACT_VALID,
    CONTRACT_MASS_MISMATCH,
    CONTRACT_POSITION_OFFSET,
    CONTRACT_VELOCITY_OFFSET
};

static int fail(const char* message) {
    fprintf(stderr, "%s\n", message);
    return EXIT_FAILURE;
}

static int try_registration(enum contract_case contract_case) {
    const double dt = 0.01;
    struct reb_simulation* main_sim = reb_simulation_create();
    struct reb_simulation* sub_sim = reb_simulation_create();
    if (!main_sim || !sub_sim) {
        if (main_sim) reb_simulation_free(main_sim);
        if (sub_sim) reb_simulation_free(sub_sim);
        return -2;
    }

    main_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    sub_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(main_sim) != 0 ||
        reb_bridge_set_integrator_whfast(sub_sim) != 0) {
        reb_simulation_free(main_sim);
        reb_simulation_free(sub_sim);
        return -2;
    }

    const double host_mass = contract_case == CONTRACT_MASS_MISMATCH ? 1.01 : 1.0;
    reb_simulation_add_fmt(main_sim, "m", host_mass);
    reb_simulation_add_fmt(sub_sim, "m", 0.75);
    reb_simulation_add_fmt(sub_sim, "m a e", 0.25, 0.02, 0.1);
    reb_simulation_move_to_com(sub_sim);

    if (contract_case == CONTRACT_POSITION_OFFSET) {
        for (int i = 0; i < sub_sim->N; i++) sub_sim->particles[i].x += 1.0e-3;
    }
    if (contract_case == CONTRACT_VELOCITY_OFFSET) {
        for (int i = 0; i < sub_sim->N; i++) sub_sim->particles[i].vx += 1.0e-3;
    }

    struct reb_bridge* bridge = reb_bridge_create(main_sim, dt);
    if (!bridge) {
        reb_simulation_free(main_sim);
        reb_simulation_free(sub_sim);
        return -2;
    }
    reb_bridge_set_owns_main(bridge, 1);
    const int status = reb_bridge_add_subsystem(bridge, 0, sub_sim, dt, 0);
    reb_bridge_free(bridge);
    reb_simulation_free(sub_sim);
    return status;
}

int main(void) {
    if (try_registration(CONTRACT_VALID) != 0) {
        return fail("canonical subsystem registration should accept a barycentric state");
    }
    if (try_registration(CONTRACT_MASS_MISMATCH) == 0) {
        return fail("subsystem registration accepted a host mass mismatch");
    }
    if (try_registration(CONTRACT_POSITION_OFFSET) == 0) {
        return fail("subsystem registration accepted a nonzero internal barycenter");
    }
    if (try_registration(CONTRACT_VELOCITY_OFFSET) == 0) {
        return fail("subsystem registration accepted nonzero internal COM velocity");
    }

    struct reb_simulation* main_sim = reb_simulation_create();
    struct reb_simulation* sub_sim = reb_simulation_create();
    if (!main_sim || !sub_sim) return fail("failed to create duplicate registration state");
    main_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    sub_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    reb_simulation_add_fmt(main_sim, "m", 1.0);
    reb_simulation_add_fmt(sub_sim, "m", 1.0);
    struct reb_bridge* bridge = reb_bridge_create(main_sim, 0.01);
    if (!bridge || reb_bridge_add_subsystem(bridge, 0, sub_sim, 0.01, 0) != 0 ||
        reb_bridge_add_subsystem(bridge, 0, sub_sim, 0.01, 0) == 0 ||
        reb_bridge_last_error() != REB_BRIDGE_ERROR_DUPLICATE_REGISTRATION) {
        if (bridge) reb_bridge_free(bridge); else reb_simulation_free(main_sim);
        reb_simulation_free(sub_sim);
        return fail("duplicate subsystem registration should be rejected");
    }
    reb_bridge_free(bridge);
    reb_simulation_free(sub_sim);
    return EXIT_SUCCESS;
}
