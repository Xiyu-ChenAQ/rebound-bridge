/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "rebound_bridge_examples.h"

#include <stdio.h>

static void examples_error(const char* msg) {
    fprintf(stderr, "[rebound_bridge examples] %s\n", msg);
}

struct reb_simulation* reb_bridge_make_main_sun_emb(void) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    reb_simulation_add_fmt(sim, "m", 1.0);
    reb_simulation_add_fmt(sim, "m a e", 3.0e-6, 1.0, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

struct reb_simulation* reb_bridge_make_sub_earth_moon(void) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    const double m_earth = 3.0e-6 * 0.987;
    const double m_moon = 3.0e-6 * 0.013;
    reb_simulation_add_fmt(sim, "m", m_earth);
    reb_simulation_add_fmt(sim, "m a e", m_moon, 0.00257, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

struct reb_bridge* reb_bridge_create_earth_moon(double dt_outer, int sub_ratio) {
    if (sub_ratio <= 0) {
        examples_error("sub_ratio must be positive");
        return NULL;
    }

    struct reb_simulation* main_sim = reb_bridge_make_main_sun_emb();
    struct reb_simulation* sub_sim = reb_bridge_make_sub_earth_moon();
    if (!main_sim || !sub_sim) {
        if (main_sim) reb_simulation_free(main_sim);
        if (sub_sim) reb_simulation_free(sub_sim);
        return NULL;
    }

    struct reb_bridge* bridge = reb_bridge_create(main_sim, dt_outer);
    if (!bridge) {
        reb_simulation_free(main_sim);
        reb_simulation_free(sub_sim);
        return NULL;
    }

    reb_bridge_set_owns_main(bridge, 1);
    if (reb_bridge_add_subsystem(bridge, 1, sub_sim, dt_outer / (double)sub_ratio, 1) != 0) {
        reb_bridge_free(bridge);
        return NULL;
    }
    return bridge;
}
