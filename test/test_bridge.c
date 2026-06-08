/* SPDX-License-Identifier: GPL-3.0-only */

/**
 * test_bridge.c - Unit tests for rebound_bridge.
 *
 * REBOUND-style C unit tests: each check is an assert().  The program
 * aborts (nonzero exit) on the first failed invariant and exits 0 only
 * if every test passes.  Wired into ctest via add_test() in CMake.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "rebound_bridge.h"

/* Relative energy drift of a single reb_simulation. */
static double sim_energy(struct reb_simulation* sim) {
    double e_kin = 0.0;
    double e_pot = 0.0;
    for (unsigned int i = 0; i < sim->N; i++) {
        const struct reb_particle pi = sim->particles[i];
        e_kin += 0.5 * pi.m * (pi.vx * pi.vx + pi.vy * pi.vy + pi.vz * pi.vz);
        for (unsigned int j = i + 1; j < sim->N; j++) {
            const struct reb_particle pj = sim->particles[j];
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            e_pot -= sim->G * pi.m * pj.m / sqrt(dx * dx + dy * dy + dz * dz);
        }
    }
    return e_kin + e_pot;
}

/* A bridge with no subsystems must reduce to plain main-sim evolution. */
static void test_create_and_free(void) {
    struct reb_simulation* main_sim = reb_bridge_make_main_sun_emb();
    assert(main_sim != NULL);

    struct reb_bridge* bridge = reb_bridge_create(main_sim, 1.0e-3);
    assert(bridge != NULL);
    reb_bridge_set_owns_main(bridge, 1);

    assert(reb_bridge_subsystem_count(bridge) == 0);
    assert(reb_bridge_main_sim(bridge) == main_sim);

    reb_bridge_free(bridge);
    printf("  [pass] test_create_and_free\n");
}

/* add_subsystem must reject malformed arguments and accept a valid one. */
static void test_add_subsystem_validation(void) {
    struct reb_simulation* main_sim = reb_bridge_make_main_sun_emb();
    assert(main_sim != NULL);
    const double dt_outer = 1.0e-3;
    struct reb_bridge* bridge = reb_bridge_create(main_sim, dt_outer);
    assert(bridge != NULL);
    reb_bridge_set_owns_main(bridge, 1);

    struct reb_simulation* sub = reb_bridge_make_sub_earth_moon();
    assert(sub != NULL);

    /* host_index out of bounds -> error, no subsystem added. */
    assert(reb_bridge_add_subsystem(bridge, 9999, sub, dt_outer / 10.0, 0) != 0);
    assert(reb_bridge_subsystem_count(bridge) == 0);

    /* dt_inner > dt_outer -> error. */
    assert(reb_bridge_add_subsystem(bridge, 1, sub, dt_outer * 2.0, 0) != 0);

    /* non-integer dt_outer/dt_inner ratio -> error. */
    assert(reb_bridge_add_subsystem(bridge, 1, sub, dt_outer / 3.3, 0) != 0);

    /* valid registration (bridge takes ownership of sub). */
    assert(reb_bridge_add_subsystem(bridge, 1, sub, dt_outer / 10.0, 1) == 0);
    assert(reb_bridge_subsystem_count(bridge) == 1);

    reb_bridge_free(bridge);
    printf("  [pass] test_add_subsystem_validation\n");
}

/* Symplectic invariant: energy drift stays bounded and small over many steps. */
static void test_energy_conservation(void) {
    struct reb_bridge* bridge = reb_bridge_create_earth_moon(1.0e-3, 20);
    assert(bridge != NULL);

    struct reb_simulation* main_sim = reb_bridge_main_sim(bridge);
    const double e0 = sim_energy(main_sim);

    double max_rel = 0.0;
    for (int i = 0; i < 2000; i++) {
        assert(reb_bridge_advance(bridge, 1.0e-3) == 0);
        const double rel = fabs((sim_energy(main_sim) - e0) / e0);
        if (rel > max_rel) max_rel = rel;
    }

    /* 2nd-order symplectic KDK: bounded, well under 1e-6 for this config. */
    assert(max_rel < 1.0e-6);

    reb_bridge_free(bridge);
    printf("  [pass] test_energy_conservation (max rel drift = %.3e)\n", max_rel);
}

int main(void) {
    printf("rebound_bridge unit tests\n");
    test_create_and_free();
    test_add_subsystem_validation();
    test_energy_conservation();
    printf("ALL TESTS PASSED\n");
    return 0;
}
