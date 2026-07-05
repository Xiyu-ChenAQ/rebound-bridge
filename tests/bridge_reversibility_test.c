/* SPDX-License-Identifier: GPL-3.0-or-later */

#include "rebound_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct particle_state {
    double x, y, z;
    double vx, vy, vz;
};

static int fail(const char* msg) {
    fprintf(stderr, "%s\n", msg);
    return EXIT_FAILURE;
}

static struct particle_state save_particle(const struct reb_particle* p) {
    struct particle_state s;
    s.x = p->x;
    s.y = p->y;
    s.z = p->z;
    s.vx = p->vx;
    s.vy = p->vy;
    s.vz = p->vz;
    return s;
}

static double particle_delta(const struct particle_state* a, const struct reb_particle* b) {
    double max_delta = 0.0;
    const double values_a[] = {a->x, a->y, a->z, a->vx, a->vy, a->vz};
    const double values_b[] = {b->x, b->y, b->z, b->vx, b->vy, b->vz};
    for (int i = 0; i < 6; i++) {
        const double delta = fabs(values_a[i] - values_b[i]);
        if (delta > max_delta) max_delta = delta;
    }
    return max_delta;
}

static struct reb_simulation* make_main_sim(double dt) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    reb_simulation_add_fmt(sim, "m", 1.0);
    reb_simulation_add_fmt(sim, "m a e", 3.0e-6, 1.0, 0.0167);
    reb_simulation_add_fmt(sim, "m a e", 9.5e-4, 5.2, 0.048);
    reb_simulation_move_to_com(sim);
    sim->dt = dt;
    return sim;
}

static struct reb_simulation* make_sub_sim(double dt) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    const double m_earth = 3.0e-6 * 0.987;
    const double m_moon = 3.0e-6 * 0.013;
    reb_simulation_add_fmt(sim, "m", m_earth);
    reb_simulation_add_fmt(sim, "m a e", m_moon, 0.00257, 0.055);
    reb_simulation_move_to_com(sim);
    sim->dt = dt;
    return sim;
}

int main(void) {
    const double dt_outer = 1.0 / 365.25;
    const double dt_inner = dt_outer / 4.0;
    const int steps = 8;

    struct reb_simulation* main_sim = make_main_sim(dt_outer);
    struct reb_simulation* sub_sim = make_sub_sim(dt_inner);
    if (!main_sim || !sub_sim) {
        if (main_sim) reb_simulation_free(main_sim);
        if (sub_sim) reb_simulation_free(sub_sim);
        return fail("failed to create simulations");
    }

    struct reb_bridge* bridge = reb_bridge_create(main_sim, dt_outer);
    if (!bridge) {
        reb_simulation_free(main_sim);
        reb_simulation_free(sub_sim);
        return fail("failed to create bridge");
    }
    reb_bridge_set_owns_main(bridge, 1);
    if (reb_bridge_add_subsystem(bridge, 1, sub_sim, dt_inner, 1) != 0) {
        reb_bridge_free(bridge);
        return fail("failed to add subsystem");
    }

    struct particle_state main_initial[3];
    struct particle_state sub_initial[2];
    for (int i = 0; i < main_sim->N; i++) {
        main_initial[i] = save_particle(&main_sim->particles[i]);
    }
    for (int i = 0; i < sub_sim->N; i++) {
        sub_initial[i] = save_particle(&sub_sim->particles[i]);
    }
    const double t0 = main_sim->t;

    for (int i = 0; i < steps; i++) {
        if (reb_bridge_step(bridge, dt_outer) != 0) {
            reb_bridge_free(bridge);
            return fail("forward bridge step failed");
        }
    }
    for (int i = 0; i < steps; i++) {
        if (reb_bridge_step(bridge, -dt_outer) != 0) {
            reb_bridge_free(bridge);
            return fail("backward bridge step failed");
        }
    }

    double max_delta = fabs(main_sim->t - t0);
    for (int i = 0; i < main_sim->N; i++) {
        const double delta = particle_delta(&main_initial[i], &main_sim->particles[i]);
        if (delta > max_delta) max_delta = delta;
    }
    for (int i = 0; i < sub_sim->N; i++) {
        const double delta = particle_delta(&sub_initial[i], &sub_sim->particles[i]);
        if (delta > max_delta) max_delta = delta;
    }

    reb_bridge_free(bridge);

    if (max_delta > 1.0e-10) {
        fprintf(stderr, "forward-backward state drift too large: %.17g\n", max_delta);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
