/* SPDX-License-Identifier: GPL-3.0-or-later */

#include "rebound_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int fail(const char* msg) {
    fprintf(stderr, "%s\n", msg);
    return EXIT_FAILURE;
}

static struct reb_simulation* make_whfast_sim(double dt) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    sim->dt = dt;
    return sim;
}

static void add_cartesian(
    struct reb_simulation* sim,
    double mass,
    double x,
    double y,
    double vx,
    double vy
) {
    struct reb_particle p = {0};
    p.m = mass;
    p.x = x;
    p.y = y;
    p.vx = vx;
    p.vy = vy;
    reb_simulation_add(sim, p);
}

static double simulation_state_delta(
    const struct reb_simulation* a,
    const struct reb_simulation* b
) {
    if (!a || !b || a->N != b->N) return INFINITY;

    double max_delta = fabs(a->t - b->t);
    for (int i = 0; i < a->N; i++) {
        const struct reb_particle* pa = &a->particles[i];
        const struct reb_particle* pb = &b->particles[i];
        const double values_a[] = {pa->x, pa->y, pa->z, pa->vx, pa->vy, pa->vz};
        const double values_b[] = {pb->x, pb->y, pb->z, pb->vx, pb->vy, pb->vz};
        for (int j = 0; j < 6; j++) {
            const double delta = fabs(values_a[j] - values_b[j]);
            if (delta > max_delta) max_delta = delta;
        }
    }
    return max_delta;
}

static int test_no_subsystem_matches_main_integrator(void) {
    const double dt = 1.0 / 365.25;
    const double duration = 32.0 * dt;
    struct reb_simulation* main_sim = make_whfast_sim(dt);
    if (!main_sim) return fail("failed to create no-subsystem main simulation");

    reb_simulation_add_fmt(main_sim, "m", 1.0);
    reb_simulation_add_fmt(main_sim, "m a e", 9.5e-4, 5.2, 0.048);
    reb_simulation_move_to_com(main_sim);

    struct reb_simulation* direct = reb_simulation_copy(main_sim);
    if (!direct) {
        reb_simulation_free(main_sim);
        return fail("failed to copy no-subsystem reference simulation");
    }

    struct reb_bridge* bridge = reb_bridge_create(main_sim, dt);
    if (!bridge) {
        reb_simulation_free(main_sim);
        reb_simulation_free(direct);
        return fail("failed to create no-subsystem bridge");
    }
    reb_bridge_set_owns_main(bridge, 1);

    const int bridge_status = reb_bridge_advance(bridge, duration);
    const enum REB_STATUS direct_status = reb_simulation_integrate(direct, duration);
    const double delta = simulation_state_delta(main_sim, direct);
    const double bridge_time = main_sim->t;
    const double direct_time = direct->t;

    reb_bridge_free(bridge);
    reb_simulation_free(direct);

    if (bridge_status != 0 || direct_status != REB_STATUS_SUCCESS) {
        return fail("no-subsystem integrations should succeed");
    }
    if (delta > 1.0e-14) {
        fprintf(
            stderr,
            "no-subsystem state delta too large: %.17g (bridge_t=%.17g direct_t=%.17g)\n",
            delta,
            bridge_time,
            direct_time
        );
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

static int test_no_external_source_matches_independent_integrators(void) {
    const double dt_outer = 1.0 / 365.25;
    const double dt_inner = dt_outer / 8.0;
    const int outer_steps = 16;

    struct reb_simulation* main_sim = make_whfast_sim(dt_outer);
    struct reb_simulation* sub_sim = make_whfast_sim(dt_inner);
    if (!main_sim || !sub_sim) {
        if (main_sim) reb_simulation_free(main_sim);
        if (sub_sim) reb_simulation_free(sub_sim);
        return fail("failed to create isolated subsystem simulations");
    }

    add_cartesian(main_sim, 1.0, 0.0, 0.0, 0.0, 0.0);
    reb_simulation_add_fmt(sub_sim, "m", 0.75);
    reb_simulation_add_fmt(sub_sim, "m a e", 0.25, 0.02, 0.1);
    reb_simulation_move_to_com(sub_sim);

    struct reb_simulation* direct_main = reb_simulation_copy(main_sim);
    struct reb_simulation* direct_sub = reb_simulation_copy(sub_sim);
    if (!direct_main || !direct_sub) {
        if (direct_main) reb_simulation_free(direct_main);
        if (direct_sub) reb_simulation_free(direct_sub);
        reb_simulation_free(main_sim);
        reb_simulation_free(sub_sim);
        return fail("failed to copy isolated subsystem references");
    }

    struct reb_bridge* bridge = reb_bridge_create(main_sim, dt_outer);
    if (!bridge) {
        reb_simulation_free(main_sim);
        reb_simulation_free(sub_sim);
        reb_simulation_free(direct_main);
        reb_simulation_free(direct_sub);
        return fail("failed to create isolated subsystem bridge");
    }
    reb_bridge_set_owns_main(bridge, 1);
    if (reb_bridge_add_subsystem(bridge, 0, sub_sim, dt_inner, 1) != 0) {
        reb_bridge_free(bridge);
        reb_simulation_free(direct_main);
        reb_simulation_free(direct_sub);
        return fail("failed to register isolated subsystem");
    }

    int bridge_status = 0;
    enum REB_STATUS main_status = REB_STATUS_SUCCESS;
    enum REB_STATUS sub_status = REB_STATUS_SUCCESS;
    for (int i = 0; i < outer_steps; i++) {
        const double target = main_sim->t + dt_outer;
        reb_simulation_synchronize(direct_main);
        reb_simulation_synchronize(direct_sub);
        bridge_status = reb_bridge_step(bridge, dt_outer);
        main_status = reb_simulation_integrate(direct_main, target);
        sub_status = reb_simulation_integrate(direct_sub, target);
        direct_sub->t = direct_main->t;
        reb_simulation_synchronize(direct_main);
        reb_simulation_synchronize(direct_sub);
        if (bridge_status != 0 ||
            main_status != REB_STATUS_SUCCESS ||
            sub_status != REB_STATUS_SUCCESS) {
            break;
        }
    }
    const double main_delta = simulation_state_delta(main_sim, direct_main);
    const double sub_delta = simulation_state_delta(sub_sim, direct_sub);

    reb_bridge_free(bridge);
    reb_simulation_free(direct_main);
    reb_simulation_free(direct_sub);

    if (bridge_status != 0 ||
        main_status != REB_STATUS_SUCCESS ||
        sub_status != REB_STATUS_SUCCESS) {
        return fail("isolated subsystem integrations should succeed");
    }
    if (main_delta > 1.0e-14 || sub_delta > 1.0e-14) {
        fprintf(
            stderr,
            "isolated subsystem state delta too large: main=%.17g sub=%.17g\n",
            main_delta,
            sub_delta
        );
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

static int measure_subsystem_kick(double source_mass, double* kick_out) {
    const double dt = 1.0e-3;
    struct reb_simulation* main_sim = make_whfast_sim(dt);
    struct reb_simulation* sub_sim = make_whfast_sim(dt);
    if (!main_sim || !sub_sim) {
        if (main_sim) reb_simulation_free(main_sim);
        if (sub_sim) reb_simulation_free(sub_sim);
        return -1;
    }

    add_cartesian(main_sim, 1.0, 0.0, 0.0, 0.0, 0.0);
    add_cartesian(main_sim, source_mass, 10.0, 1.0, 0.0, 0.0);
    add_cartesian(sub_sim, 0.75, -0.01, 0.0, 0.0, 0.0);
    add_cartesian(sub_sim, 0.25, 0.03, 0.0, 0.0, 0.0);

    struct reb_bridge* bridge = reb_bridge_create(main_sim, dt);
    if (!bridge) {
        reb_simulation_free(main_sim);
        reb_simulation_free(sub_sim);
        return -1;
    }
    reb_bridge_set_owns_main(bridge, 1);
    if (reb_bridge_add_subsystem(bridge, 0, sub_sim, dt, 1) != 0) {
        reb_bridge_free(bridge);
        return -1;
    }

    const int status = reb_bridge_apply_cross_kick(bridge, dt);
    if (status == 0) {
        const struct reb_particle* p = &sub_sim->particles[0];
        *kick_out = sqrt(p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
    }
    reb_bridge_free(bridge);
    return status;
}

static int test_cross_kick_vanishes_with_coupling_strength(void) {
    double zero_kick = 0.0;
    double strong_kick = 0.0;
    double weak_kick = 0.0;
    if (measure_subsystem_kick(0.0, &zero_kick) != 0 ||
        measure_subsystem_kick(1.0, &strong_kick) != 0 ||
        measure_subsystem_kick(1.0e-3, &weak_kick) != 0) {
        return fail("failed to measure limiting cross-kicks");
    }
    if (zero_kick != 0.0) {
        fprintf(stderr, "zero-mass source produced a nonzero kick: %.17g\n", zero_kick);
        return EXIT_FAILURE;
    }
    if (strong_kick <= 0.0 || weak_kick <= 0.0) {
        return fail("nonzero external sources should produce nonzero tidal kicks");
    }

    const double ratio = weak_kick / strong_kick;
    if (fabs(ratio - 1.0e-3) > 1.0e-12) {
        fprintf(stderr, "cross-kick did not scale linearly with source mass: %.17g\n", ratio);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int main(void) {
    if (test_no_subsystem_matches_main_integrator() != EXIT_SUCCESS) return EXIT_FAILURE;
    if (test_no_external_source_matches_independent_integrators() != EXIT_SUCCESS) return EXIT_FAILURE;
    if (test_cross_kick_vanishes_with_coupling_strength() != EXIT_SUCCESS) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
