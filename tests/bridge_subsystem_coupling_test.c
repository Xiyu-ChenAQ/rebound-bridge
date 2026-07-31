/* SPDX-License-Identifier: GPL-3.0-or-later */

#include "rebound_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int fail(const char* msg) {
    fprintf(stderr, "%s\n", msg);
    return EXIT_FAILURE;
}

static void add_particle(struct reb_simulation* sim, double m, double x, double y, double z) {
    struct reb_particle p = {0};
    p.m = m;
    p.x = x;
    p.y = y;
    p.z = z;
    reb_simulation_add(sim, p);
}

static struct reb_simulation* make_subsystem_a(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    add_particle(sim, 2.0, -0.01, 0.00, 0.0);
    add_particle(sim, 1.0,  0.02, 0.01, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_subsystem_b(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    add_particle(sim, 1.5, 0.00, -0.02, 0.0);
    add_particle(sim, 0.5, 0.03,  0.06, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_bridge* make_two_subsystem_bridge(
    enum reb_bridge_subsystem_coupling coupling,
    struct reb_simulation** sub_a_out,
    struct reb_simulation** sub_b_out
) {
    struct reb_simulation* main_sim = reb_simulation_create();
    if (!main_sim) return NULL;
    main_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(main_sim) != 0) {
        reb_simulation_free(main_sim);
        return NULL;
    }
    add_particle(main_sim, 3.0, 0.0, 0.0, 0.0);
    add_particle(main_sim, 2.0, 1.0, 0.1, 0.0);

    struct reb_bridge* bridge = reb_bridge_create(main_sim, 0.01);
    if (!bridge) {
        reb_simulation_free(main_sim);
        return NULL;
    }
    reb_bridge_set_owns_main(bridge, 1);

    struct reb_simulation* sub_a = make_subsystem_a();
    struct reb_simulation* sub_b = make_subsystem_b();
    if (!sub_a || !sub_b) {
        if (sub_a) reb_simulation_free(sub_a);
        if (sub_b) reb_simulation_free(sub_b);
        reb_bridge_free(bridge);
        return NULL;
    }

    if (reb_bridge_add_subsystem(bridge, 0, sub_a, 0.01, 1) != 0 ||
        reb_bridge_add_subsystem(bridge, 1, sub_b, 0.01, 1) != 0 ||
        reb_bridge_set_subsystem_coupling(bridge, coupling) != 0) {
        reb_bridge_free(bridge);
        return NULL;
    }

    *sub_a_out = sub_a;
    *sub_b_out = sub_b;
    return bridge;
}

static double internal_momentum_norm(const struct reb_simulation* sim) {
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    for (int i = 0; i < sim->N; i++) {
        px += sim->particles[i].m * sim->particles[i].vx;
        py += sim->particles[i].m * sim->particles[i].vy;
        pz += sim->particles[i].m * sim->particles[i].vz;
    }
    return sqrt(px * px + py * py + pz * pz);
}

static double particle_velocity_delta(const struct reb_particle* a, const struct reb_particle* b) {
    const double dx = a->vx - b->vx;
    const double dy = a->vy - b->vy;
    const double dz = a->vz - b->vz;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

int main(void) {
    if (reb_bridge_get_subsystem_coupling(NULL) != REB_BRIDGE_SUBSYSTEM_COUPLING_INVALID) {
        return fail("NULL bridge should report invalid subsystem coupling");
    }

    struct reb_simulation* mono_a = NULL;
    struct reb_simulation* mono_b = NULL;
    struct reb_bridge* mono = make_two_subsystem_bridge(
        REB_BRIDGE_SUBSYSTEM_COUPLING_MONOPOLE,
        &mono_a,
        &mono_b
    );
    if (!mono) return fail("failed to create monopole bridge");

    if (reb_bridge_get_subsystem_coupling(mono) != REB_BRIDGE_SUBSYSTEM_COUPLING_MONOPOLE) {
        reb_bridge_free(mono);
        return fail("new bridge should default to monopole subsystem coupling");
    }
    if (reb_bridge_set_subsystem_coupling(mono, (enum reb_bridge_subsystem_coupling)999) == 0) {
        reb_bridge_free(mono);
        return fail("invalid subsystem coupling should fail");
    }
    if (reb_bridge_get_subsystem_coupling(mono) != REB_BRIDGE_SUBSYSTEM_COUPLING_MONOPOLE) {
        reb_bridge_free(mono);
        return fail("invalid subsystem coupling should not change current mode");
    }

    struct reb_simulation* resolved_a = NULL;
    struct reb_simulation* resolved_b = NULL;
    struct reb_bridge* resolved = make_two_subsystem_bridge(
        REB_BRIDGE_SUBSYSTEM_COUPLING_RESOLVED,
        &resolved_a,
        &resolved_b
    );
    if (!resolved) {
        reb_bridge_free(mono);
        return fail("failed to create resolved bridge");
    }
    if (reb_bridge_get_subsystem_coupling(resolved) != REB_BRIDGE_SUBSYSTEM_COUPLING_RESOLVED) {
        reb_bridge_free(mono);
        reb_bridge_free(resolved);
        return fail("bridge should remember resolved subsystem coupling");
    }

    if (reb_bridge_apply_cross_kick(mono, 1.0e-4) != 0 ||
        reb_bridge_apply_cross_kick(resolved, 1.0e-4) != 0) {
        reb_bridge_free(mono);
        reb_bridge_free(resolved);
        return fail("cross-kicks should succeed");
    }

    const double delta_a0 = particle_velocity_delta(&mono_a->particles[0], &resolved_a->particles[0]);
    const double delta_b0 = particle_velocity_delta(&mono_b->particles[0], &resolved_b->particles[0]);
    if (delta_a0 <= 1.0e-12 || delta_b0 <= 1.0e-12) {
        reb_bridge_free(mono);
        reb_bridge_free(resolved);
        return fail("resolved coupling should change internal subsystem kicks relative to monopole coupling");
    }

    if (internal_momentum_norm(resolved_a) > 1.0e-14 ||
        internal_momentum_norm(resolved_b) > 1.0e-14) {
        reb_bridge_free(mono);
        reb_bridge_free(resolved);
        return fail("resolved coupling should preserve each subsystem internal COM momentum");
    }

    reb_bridge_free(mono);
    reb_bridge_free(resolved);
    return EXIT_SUCCESS;
}
