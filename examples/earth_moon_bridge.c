/* SPDX-License-Identifier: GPL-3.0-or-later */

/**
 * earth_moon_bridge.c - Standalone rebound_bridge validation example.
 *
 * Demonstrates convergence order and energy conservation of the bridge
 * integrator for an Earth-Moon subsystem inside a Sun-EMB main simulation.
 * Diagnostic quantities are read directly from the sub-simulation via
 * reb_bridge_sub_sim(), which is the intended public query interface.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rebound_bridge.h"
#include "rebound_bridge_examples.h"

/* Scalar distance between two particles. */
static double particle_distance(
    const struct reb_particle* a,
    const struct reb_particle* b
) {
    const double dx = a->x - b->x;
    const double dy = a->y - b->y;
    const double dz = a->z - b->z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/*
 * Two-body specific energy of the subsystem.
 * Assumes particles[0] = Earth, particles[1] = Moon.
 */
static double sub_energy(const struct reb_simulation* sub_sim) {
    const struct reb_particle* e = &sub_sim->particles[0];
    const struct reb_particle* m = &sub_sim->particles[1];
    const double dvx = e->vx - m->vx;
    const double dvy = e->vy - m->vy;
    const double dvz = e->vz - m->vz;
    const double mu  = e->m * m->m / (e->m + m->m);
    const double K   = 0.5 * mu * (dvx*dvx + dvy*dvy + dvz*dvz);
    const double r   = particle_distance(e, m);
    const double U   = -sub_sim->G * e->m * m->m / r;
    return K + U;
}

static int collect_distances(double dt_outer, double* distances, int n_samples) {
    const double t_end = 1.0;
    struct reb_bridge* bridge = reb_bridge_create_earth_moon(dt_outer, 50);
    if (!bridge) return -1;

    for (int i = 0; i < n_samples; i++) {
        const double target = t_end * (double)(i + 1) / (double)n_samples;
        struct reb_simulation* main_sim = reb_bridge_main_sim(bridge);
        if (reb_bridge_advance(bridge, target - main_sim->t) != 0) {
            reb_bridge_free(bridge);
            return -1;
        }
        struct reb_simulation* sub_sim = reb_bridge_sub_sim(bridge, 0);
        distances[i] = particle_distance(&sub_sim->particles[0], &sub_sim->particles[1]);
    }

    reb_bridge_free(bridge);
    return 0;
}

static double rms_difference(const double* a, const double* b, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        const double d = a[i] - b[i];
        sum += d * d;
    }
    return sqrt(sum / (double)n);
}

static void convergence_demo(void) {
    const int n_dt = 4;
    const int n_samples = 50;
    const double dts[] = {1.0/365.0, 1.0/730.0, 1.0/1460.0, 1.0/2920.0};
    double* results[4] = {0};

    printf("rebound_bridge - convergence test\n");
    printf("==================================\n");

    for (int i = 0; i < n_dt; i++) {
        results[i] = (double*)calloc((size_t)n_samples, sizeof(double));
        if (!results[i] || collect_distances(dts[i], results[i], n_samples) != 0) {
            printf("  dt=%g failed\n", dts[i]);
            goto cleanup;
        }
        printf("  dt=1/%-4.0f yr done\n", 1.0 / dts[i]);
    }

    printf("\n  %-16s %-14s %-8s %-8s\n", "dt", "RMS vs finest", "ratio", "order");
    printf("  -----------------------------------------------\n");
    for (int i = 0; i < n_dt - 2; i++) {
        const double e1 = rms_difference(results[i],     results[n_dt - 1], n_samples);
        const double e2 = rms_difference(results[i + 1], results[n_dt - 1], n_samples);
        const double ratio = e1 / e2;
        const double order = log(ratio) / log(2.0);
        printf("  1/%-5.0f -> 1/%-5.0f  %.3e      %.3f    %.3f\n",
               1.0 / dts[i], 1.0 / dts[i + 1], e1, ratio, order);
    }
    printf("  1/%-5.0f -> reference  %.3e\n",
           1.0 / dts[n_dt - 2],
           rms_difference(results[n_dt - 2], results[n_dt - 1], n_samples));

cleanup:
    for (int i = 0; i < n_dt; i++) free(results[i]);
}

static void diagnostics_demo(void) {
    const int n_samples = 50;
    const double dt_outer = 1.0 / 365.0;
    struct reb_bridge* bridge = reb_bridge_create_earth_moon(dt_outer, 50);
    if (!bridge) {
        printf("diagnostics failed: could not create bridge\n");
        return;
    }

    double mean_d = 0.0, mean_d2 = 0.0;
    double e0 = 0.0, max_de = 0.0;
    int have_e0 = 0;

    for (int i = 0; i < n_samples; i++) {
        const double target = (double)(i + 1) / (double)n_samples;
        struct reb_simulation* main_sim = reb_bridge_main_sim(bridge);
        if (reb_bridge_advance(bridge, target - main_sim->t) != 0) {
            printf("diagnostics failed during advance\n");
            reb_bridge_free(bridge);
            return;
        }

        struct reb_simulation* sub_sim = reb_bridge_sub_sim(bridge, 0);
        const double d = particle_distance(&sub_sim->particles[0], &sub_sim->particles[1]);
        const double E = sub_energy(sub_sim);
        if (!have_e0) { e0 = E; have_e0 = 1; }
        const double de = fabs(E - e0);
        if (de > max_de) max_de = de;
        mean_d  += d;
        mean_d2 += d * d;
    }

    mean_d  /= (double)n_samples;
    mean_d2 /= (double)n_samples;
    const double std_d = sqrt(mean_d2 - mean_d * mean_d);

    printf("\nDiagnostics (dt=1/365 yr)\n");
    printf("==========================\n");
    printf("  d_mean = %.8f AU\n", mean_d);
    printf("  d_std  = %.3e AU\n", std_d);
    printf("  dE_max = %.3e\n",    max_de);

    reb_bridge_free(bridge);
}

int main(void) {
    convergence_demo();
    diagnostics_demo();
    return 0;
}
