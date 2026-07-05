/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SOLAR_SYSTEM_FIXTURE_H
#define SOLAR_SYSTEM_FIXTURE_H

/**
 * @file  solar_system_fixture.h
 * @brief Shared solar-system initial conditions for validation programs.
 *
 * All positions are J2000 Horizons state vectors in AU/(AU/yr) (G units:
 * REB_BRIDGE_G_AU_YR_MSUN).  Masses are in solar masses.
 */

#include "rebound.h"
#include "rebound_bridge.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Body index enumerations shared by solar-system validation programs.
/* Particles in the main (barycentric) simulation when using bridge mode. */
enum ss_main_body {
    SS_MAIN_SUN              = 0,
    SS_MAIN_MERCURY          = 1,
    SS_MAIN_VENUS            = 2,
    SS_MAIN_EARTH_MOON_HOST  = 3,
    SS_MAIN_MARS             = 4,
    SS_MAIN_JUPITER_HOST     = 5,
    SS_MAIN_SATURN           = 6,
    SS_MAIN_URANUS           = 7,
    SS_MAIN_NEPTUNE          = 8
};

/* Particles in the full-resolution (WHFast/IAS15) simulation. */
enum ss_body {
    SS_BODY_SUN     = 0,
    SS_BODY_MERCURY = 1,
    SS_BODY_VENUS   = 2,
    SS_BODY_EARTH   = 3,
    SS_BODY_MOON    = 4,
    SS_BODY_MARS    = 5,
    SS_BODY_JUPITER = 6,
    SS_BODY_IO      = 7,
    SS_BODY_EUROPA  = 8,
    SS_BODY_GANYMEDE= 9,
    SS_BODY_SATURN  = 10,
    SS_BODY_URANUS  = 11,
    SS_BODY_NEPTUNE = 12,
    SS_BODY_COUNT   = 13
};

enum ss_em_body  { SS_EM_EARTH = 0, SS_EM_MOON = 1 };
enum ss_jov_body {
    SS_JOV_JUPITER = 0,
    SS_JOV_IO      = 1,
    SS_JOV_EUROPA  = 2,
    SS_JOV_GANYMEDE= 3
};

/*
 * Fill bodies[SS_BODY_COUNT] with J2000 Horizons state vectors.
 * The array must be pre-allocated with at least SS_BODY_COUNT elements.
 */
void ss_populate_bodies(struct reb_particle bodies[SS_BODY_COUNT]);

/* Add all 13 bodies to sim (solar system, full resolution). */
void ss_add_all_bodies(struct reb_simulation* sim);

/*
 * Create a full solar-system WHFast simulation (all 13 bodies).
 * dt_days sets the timestep.  Caller owns the returned simulation.
 */
struct reb_simulation* ss_make_whfast(double dt_days);

/*
 * Create a full solar-system IAS15 simulation (all 13 bodies).
 * Caller owns the returned simulation.
 */
struct reb_simulation* ss_make_ias15(void);

#ifdef __cplusplus
}
#endif

#endif /* SOLAR_SYSTEM_FIXTURE_H */
