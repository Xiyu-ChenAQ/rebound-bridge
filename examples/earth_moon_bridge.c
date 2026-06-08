/* SPDX-License-Identifier: GPL-3.0-only */

/**
 * earth_moon_bridge.c - Standalone rebound_bridge validation example.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rebound_bridge.h"

static double distance3(
    const struct reb_bridge_body_state* a,
    const struct reb_bridge_body_state* b
) {
    const double dx = a->x - b->x;
    const double dy = a->y - b->y;
    const double dz = a->z - b->z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

static double sub_energy(const struct reb_bridge_snapshot* s) {
    const double me = s->member0.m;
    const double mm = s->member1.m;
    const double dvx = s->member0.vx - s->member1.vx;
    const double dvy = s->member0.vy - s->member1.vy;
    const double dvz = s->member0.vz - s->member1.vz;
    const double mu = me * mm / (me + mm);
    const double K = 0.5 * mu * (dvx*dvx + dvy*dvy + dvz*dvz);
    const double U = -REB_BRIDGE_G_AU_YR_MSUN * me * mm / s->subsystem_internal_distance;
    return K + U;
}

static int collect_distances(double dt_outer, double* distances, int n_samples) {
    const double t_end = 1.0;
    struct reb_bridge* bridge = reb_bridge_create_earth_moon(dt_outer, 50);
    if (!bridge) return -1;

    for (int i = 0; i < n_samples; i++) {
        const double target = t_end * (double)(i + 1) / (double)n_samples;
        struct reb_simulation* main_sim = reb_bridge_main_sim(bridge);
        const double duration = target - main_sim->t;
        if (reb_bridge_advance(bridge, duration) != 0) {
            reb_bridge_free(bridge);
            return -1;
        }

        struct reb_bridge_snapshot snap;
        if (reb_bridge_snapshot(bridge, 0, &snap) != 0) {
            reb_bridge_free(bridge);
            return -1;
        }
        distances[i] = snap.subsystem_internal_distance;
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
        const double e1 = rms_difference(results[i], results[n_dt - 1], n_samples);
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
    for (int i = 0; i < n_dt; i++) {
        free(results[i]);
    }
}

static void diagnostics_demo(void) {
    const int n_samples = 50;
    const double dt_outer = 1.0 / 365.0;
    struct reb_bridge* bridge = reb_bridge_create_earth_moon(dt_outer, 50);
    if (!bridge) {
        printf("diagnostics failed: could not create bridge\n");
        return;
    }

    double mean_d = 0.0;
    double mean_d2 = 0.0;
    double e0 = 0.0;
    double max_de = 0.0;
    int have_e0 = 0;

    for (int i = 0; i < n_samples; i++) {
        const double target = (double)(i + 1) / (double)n_samples;
        struct reb_simulation* main_sim = reb_bridge_main_sim(bridge);
        if (reb_bridge_advance(bridge, target - main_sim->t) != 0) {
            printf("diagnostics failed during advance\n");
            reb_bridge_free(bridge);
            return;
        }

        struct reb_bridge_snapshot s;
        if (reb_bridge_snapshot(bridge, 0, &s) != 0) {
            printf("diagnostics failed during snapshot\n");
            reb_bridge_free(bridge);
            return;
        }

        const double d = distance3(&s.member0, &s.member1);
        const double E = sub_energy(&s);
        if (!have_e0) {
            e0 = E;
            have_e0 = 1;
        }
        const double de = fabs(E - e0);
        if (de > max_de) max_de = de;
        mean_d += d;
        mean_d2 += d * d;
    }

    mean_d /= (double)n_samples;
    mean_d2 /= (double)n_samples;
    const double std_d = sqrt(mean_d2 - mean_d * mean_d);

    printf("\nDiagnostics (dt=1/365 yr)\n");
    printf("==========================\n");
    printf("  d_mean = %.8f AU\n", mean_d);
    printf("  d_std  = %.3e AU\n", std_d);
    printf("  dE_max = %.3e\n", max_de);

    reb_bridge_free(bridge);
}

int main(void) {
    convergence_demo();
    diagnostics_demo();
    return 0;
}
