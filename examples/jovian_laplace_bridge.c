/* SPDX-License-Identifier: GPL-3.0-or-later */

/**
 * jovian_laplace_bridge.c - Standalone rebound_bridge Laplace resonance example.
 *
 * This example keeps the Jovian system resolved in a bridge subsystem while the
 * main simulation carries the Sun and a single Jupiter-system host particle.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rebound_bridge.h"

enum jovian_body {
    BODY_JUPITER = 0,
    BODY_IO = 1,
    BODY_EUROPA = 2,
    BODY_GANYMEDE = 3
};

static double wrap_degrees(double angle) {
    angle = fmod(angle, 360.0);
    if (angle < 0.0) angle += 360.0;
    return angle;
}

static struct reb_simulation* make_main_sun_jupiter(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    const double m_jupiter_system = 9.547919e-4 + 2.388e-8 + 2.527e-8 + 7.804e-8;
    reb_simulation_add_fmt(sim, "m", 1.0);
    reb_simulation_add_fmt(sim, "m a e", m_jupiter_system, 5.2044, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

static void add_moon(
    struct reb_simulation* sim,
    double mass,
    double x,
    double y,
    double z,
    double vx,
    double vy,
    double vz
) {
    reb_simulation_add_fmt(sim, "m x y z vx vy vz", mass, x, y, z, vx, vy, vz);
}

static struct reb_simulation* make_sub_jovian_laplace(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    const double m_jupiter = 9.547919e-4;
    const double m_io = 2.388e-8;
    const double m_europa = 2.527e-8;
    const double m_ganymede = 7.804e-8;

    /*
     * Demo-tuned initial state from JPL Horizons for 2026-Jan-01 00:00 TDB,
     * center 500@599 (Jupiter), geometric vectors, AU/day converted to AU/yr.
     * This is fixed data for a stable example, not a general ephemeris loader.
     */
    reb_simulation_add_fmt(sim, "m", m_jupiter);
    add_moon(
        sim,
        m_io,
        2.48120485765584205e-03,
        -1.34556438159649208e-03,
        -1.31630007347439604e-05,
        1.75341752994026812e+00,
        3.20163492363371871e+00,
        1.39542430489110164e-01
    );
    add_moon(
        sim,
        m_europa,
        5.51813635549570501e-04,
        -4.47551321186325738e-03,
        -1.27610904592036212e-04,
        2.85755306406471643e+00,
        3.72439324002732552e-01,
        7.16871253077594472e-02
    );
    add_moon(
        sim,
        m_ganymede,
        6.76397567656043287e-03,
        -2.29852901298470890e-03,
        1.01873701583559195e-05,
        7.40898140660848692e-01,
        2.17447712709783314e+00,
        9.37799108604124437e-02
    );
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_bridge* create_jovian_laplace_bridge(double dt_outer, int sub_ratio) {
    if (sub_ratio <= 0) return NULL;

    struct reb_simulation* main_sim = make_main_sun_jupiter();
    struct reb_simulation* sub_sim = make_sub_jovian_laplace();
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

static struct reb_orbit moon_orbit(const struct reb_simulation* sub_sim, int moon_index) {
    return reb_orbit_from_particle(
        sub_sim->G,
        sub_sim->particles[moon_index],
        sub_sim->particles[BODY_JUPITER]
    );
}

static double laplace_angle_degrees(const struct reb_simulation* sub_sim) {
    const struct reb_orbit io = moon_orbit(sub_sim, BODY_IO);
    const struct reb_orbit europa = moon_orbit(sub_sim, BODY_EUROPA);
    const struct reb_orbit ganymede = moon_orbit(sub_sim, BODY_GANYMEDE);
    return wrap_degrees((io.l - 3.0 * europa.l + 2.0 * ganymede.l) * 180.0 / M_PI);
}

static void print_jovian_diagnostics(const char* label, const struct reb_simulation* sub_sim) {
    const struct reb_orbit io = moon_orbit(sub_sim, BODY_IO);
    const struct reb_orbit europa = moon_orbit(sub_sim, BODY_EUROPA);
    const struct reb_orbit ganymede = moon_orbit(sub_sim, BODY_GANYMEDE);

    printf("%s\n", label);
    printf("  a_io       = %.8f AU\n", io.a);
    printf("  a_europa   = %.8f AU\n", europa.a);
    printf("  a_ganymede = %.8f AU\n", ganymede.a);
    printf("  P_europa / P_io       = %.8f\n", europa.P / io.P);
    printf("  P_ganymede / P_europa = %.8f\n", ganymede.P / europa.P);
    printf("  phi_L = lambda_io - 3 lambda_europa + 2 lambda_ganymede = %.6f deg\n",
           laplace_angle_degrees(sub_sim));
}

int main(void) {
    const double dt_outer = 1.0 / 3650.0;
    const int sub_ratio = 20;
    const double duration = 30.0 / 365.25;

    struct reb_bridge* bridge = create_jovian_laplace_bridge(dt_outer, sub_ratio);
    if (!bridge) {
        fprintf(stderr, "failed to create Jovian Laplace bridge\n");
        return EXIT_FAILURE;
    }

    struct reb_simulation* sub_sim = reb_bridge_sub_sim(bridge, 0);
    printf("rebound_bridge - Jovian Laplace resonance example\n");
    printf("=================================================\n");
    print_jovian_diagnostics("Initial state", sub_sim);

    if (reb_bridge_advance(bridge, duration) != 0) {
        reb_bridge_free(bridge);
        return EXIT_FAILURE;
    }

    printf("\nAfter %.2f days\n", duration * 365.25);
    print_jovian_diagnostics("Final state", sub_sim);

    reb_bridge_free(bridge);
    return EXIT_SUCCESS;
}
