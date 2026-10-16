/* SPDX-License-Identifier: GPL-3.0-or-later */

#include "rebound_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double DT_VALUES_DAYS[] = {0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
static const double COUPLING_SCALES[] = {0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 5.0};
static const double DURATION_DAYS = 128.0;
static const double REFERENCE_DT_DAYS = 0.03125;
static const double FINE_REFERENCE_DT_DAYS = 0.015625;
static const double M_EARTH = 3.0e-6;
static const double M_MOON = 3.7e-8;
static const double MOON_A = 0.00257;

enum split_method {
    SPLIT_ORDER2,
    SPLIT_LF42
};

struct experiment_state {
    struct reb_bridge* bridge;
    struct reb_simulation* main_sim;
    struct reb_simulation* sub_sim;
};

static struct experiment_state make_state(void) {
    struct experiment_state state = {0};
    const double subsystem_mass = M_EARTH + M_MOON;

    state.main_sim = reb_simulation_create();
    state.sub_sim = reb_simulation_create();
    if (!state.main_sim || !state.sub_sim) return state;

    state.main_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    state.sub_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_ias15(state.main_sim) != 0 ||
        reb_bridge_set_integrator_ias15(state.sub_sim) != 0) {
        return state;
    }
    state.main_sim->ri_ias15.epsilon = 1.0e-13;
    state.sub_sim->ri_ias15.epsilon = 1.0e-13;
    state.main_sim->exact_finish_time = 1;
    state.sub_sim->exact_finish_time = 1;

    reb_simulation_add_fmt(state.main_sim, "m", 1.0);
    reb_simulation_add_fmt(state.main_sim, "m a e", subsystem_mass, 1.0, 0.02);
    reb_simulation_move_to_com(state.main_sim);

    reb_simulation_add_fmt(state.sub_sim, "m", M_EARTH);
    reb_simulation_add_fmt(state.sub_sim, "m a e", M_MOON, MOON_A, 0.05);
    reb_simulation_move_to_com(state.sub_sim);

    state.bridge = reb_bridge_create(state.main_sim, 1.0 / 365.25);
    if (!state.bridge) return state;
    reb_bridge_set_owns_main(state.bridge, 1);
    if (reb_bridge_add_subsystem(state.bridge, 1, state.sub_sim, 1.0 / 365.25, 1) != 0) {
        return state;
    }
    return state;
}

static void free_state(struct experiment_state* state) {
    if (state->bridge) {
        reb_bridge_free(state->bridge);
    } else {
        if (state->main_sim) reb_simulation_free(state->main_sim);
        if (state->sub_sim) reb_simulation_free(state->sub_sim);
    }
    *state = (struct experiment_state){0};
}

static int integrate_a(struct experiment_state* state, double dt) {
    const double target = state->main_sim->t + dt;
    if (reb_simulation_integrate(state->sub_sim, target) != REB_STATUS_SUCCESS) return -1;
    if (reb_simulation_integrate(state->main_sim, target) != REB_STATUS_SUCCESS) return -1;
    return 0;
}

static int scaled_kick(struct experiment_state* state, double dt, double coupling_scale) {
    if (state->main_sim->N != 2 || state->sub_sim->N != 2) return -1;
    double main_velocity[2][3];
    double sub_velocity[2][3];
    for (int i = 0; i < 2; i++) {
        const struct reb_particle* main_particle = &state->main_sim->particles[i];
        const struct reb_particle* sub_particle = &state->sub_sim->particles[i];
        main_velocity[i][0] = main_particle->vx;
        main_velocity[i][1] = main_particle->vy;
        main_velocity[i][2] = main_particle->vz;
        sub_velocity[i][0] = sub_particle->vx;
        sub_velocity[i][1] = sub_particle->vy;
        sub_velocity[i][2] = sub_particle->vz;
    }

    if (reb_bridge_apply_cross_kick(state->bridge, dt) != 0) return -1;

    for (int i = 0; i < 2; i++) {
        struct reb_particle* main_particle = &state->main_sim->particles[i];
        struct reb_particle* sub_particle = &state->sub_sim->particles[i];
        main_particle->vx = main_velocity[i][0] + coupling_scale * (main_particle->vx - main_velocity[i][0]);
        main_particle->vy = main_velocity[i][1] + coupling_scale * (main_particle->vy - main_velocity[i][1]);
        main_particle->vz = main_velocity[i][2] + coupling_scale * (main_particle->vz - main_velocity[i][2]);
        sub_particle->vx = sub_velocity[i][0] + coupling_scale * (sub_particle->vx - sub_velocity[i][0]);
        sub_particle->vy = sub_velocity[i][1] + coupling_scale * (sub_particle->vy - sub_velocity[i][1]);
        sub_particle->vz = sub_velocity[i][2] + coupling_scale * (sub_particle->vz - sub_velocity[i][2]);
    }
    return 0;
}

static int split_step(
    struct experiment_state* state,
    enum split_method method,
    double dt,
    double coupling_scale
) {
    if (method == SPLIT_ORDER2) {
        if (scaled_kick(state, 0.5 * dt, coupling_scale) != 0) return -1;
        if (integrate_a(state, dt) != 0) return -1;
        return scaled_kick(state, 0.5 * dt, coupling_scale);
    }

    const double a = 0.211324865405187117745425609749;
    if (integrate_a(state, a * dt) != 0) return -1;
    if (scaled_kick(state, 0.5 * dt, coupling_scale) != 0) return -1;
    if (integrate_a(state, (1.0 - 2.0 * a) * dt) != 0) return -1;
    if (scaled_kick(state, 0.5 * dt, coupling_scale) != 0) return -1;
    return integrate_a(state, a * dt);
}

static int run_split(
    struct experiment_state* state,
    enum split_method method,
    double dt_days,
    double coupling_scale
) {
    const double dt = dt_days / 365.25;
    const int steps = (int)llround(DURATION_DAYS / dt_days);
    if (fabs((double)steps * dt_days - DURATION_DAYS) > 1.0e-12) return -1;
    for (int i = 0; i < steps; i++) {
        if (split_step(state, method, dt, coupling_scale) != 0) return -1;
    }
    return 0;
}

static double position_error(
    const struct experiment_state* state,
    const struct experiment_state* reference,
    int subsystem
) {
    const struct reb_simulation* sim = subsystem ? state->sub_sim : state->main_sim;
    const struct reb_simulation* ref = subsystem ? reference->sub_sim : reference->main_sim;
    double max_error = 0.0;
    for (int i = 0; i < sim->N; i++) {
        const double dx = sim->particles[i].x - ref->particles[i].x;
        const double dy = sim->particles[i].y - ref->particles[i].y;
        const double dz = sim->particles[i].z - ref->particles[i].z;
        const double error = sqrt(dx * dx + dy * dy + dz * dz);
        if (error > max_error) max_error = error;
    }
    return max_error;
}

static const char* method_name(enum split_method method) {
    return method == SPLIT_ORDER2 ? "order2" : "lf42";
}

int main(int argc, char** argv) {
    if (argc != 3 || strcmp(argv[1], "--csv") != 0) {
        fprintf(stderr, "usage: %s --csv PATH\n", argv[0]);
        return EXIT_FAILURE;
    }
    FILE* csv = fopen(argv[2], "w");
    if (!csv) {
        fprintf(stderr, "failed to open %s\n", argv[2]);
        return EXIT_FAILURE;
    }
    fprintf(
        csv,
        "coupling_scale,effective_epsilon,method,dt_days,duration_days,"
        "sub_position_error_au,main_position_error_au,max_position_error_au,reference_floor_au\n"
    );

    const double base_epsilon = (1.0 / (M_EARTH + M_MOON)) * MOON_A * MOON_A * MOON_A;
    const int scale_count = (int)(sizeof(COUPLING_SCALES) / sizeof(COUPLING_SCALES[0]));
    const int dt_count = (int)(sizeof(DT_VALUES_DAYS) / sizeof(DT_VALUES_DAYS[0]));

    for (int scale_index = 0; scale_index < scale_count; scale_index++) {
        const double coupling_scale = COUPLING_SCALES[scale_index];
        struct experiment_state reference = make_state();
        struct experiment_state fine_reference = make_state();
        if (!reference.bridge || !fine_reference.bridge ||
            run_split(&reference, SPLIT_LF42, REFERENCE_DT_DAYS, coupling_scale) != 0 ||
            run_split(&fine_reference, SPLIT_LF42, FINE_REFERENCE_DT_DAYS, coupling_scale) != 0) {
            free_state(&reference);
            free_state(&fine_reference);
            fclose(csv);
            return EXIT_FAILURE;
        }
        const double reference_floor = fmax(
            position_error(&reference, &fine_reference, 0),
            position_error(&reference, &fine_reference, 1)
        );

        for (int method_value = SPLIT_ORDER2; method_value <= SPLIT_LF42; method_value++) {
            const enum split_method method = (enum split_method)method_value;
            for (int dt_index = 0; dt_index < dt_count; dt_index++) {
                struct experiment_state state = make_state();
                if (!state.bridge ||
                    run_split(&state, method, DT_VALUES_DAYS[dt_index], coupling_scale) != 0) {
                    free_state(&state);
                    free_state(&reference);
                    free_state(&fine_reference);
                    fclose(csv);
                    return EXIT_FAILURE;
                }
                const double sub_error = position_error(&state, &fine_reference, 1);
                const double main_error = position_error(&state, &fine_reference, 0);
                fprintf(
                    csv,
                    "%.17g,%.17g,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    coupling_scale,
                    coupling_scale * base_epsilon,
                    method_name(method),
                    DT_VALUES_DAYS[dt_index],
                    DURATION_DAYS,
                    sub_error,
                    main_error,
                    fmax(sub_error, main_error),
                    reference_floor
                );
                free_state(&state);
            }
        }
        free_state(&reference);
        free_state(&fine_reference);
    }

    fclose(csv);
    return EXIT_SUCCESS;
}
