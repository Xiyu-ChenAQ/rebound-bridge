/* SPDX-License-Identifier: GPL-3.0-or-later */

#include "galilean_horizons_state.h"
#include "rebound_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum full_body {
    BODY_SUN = 0,
    BODY_JUPITER = 1,
    BODY_IO = 2,
    BODY_EUROPA = 3,
    BODY_GANYMEDE = 4,
    BODY_CALLISTO = 5,
    BODY_COUNT = 6
};

struct bridge_state {
    struct reb_simulation* main_sim;
    struct reb_simulation* sub_sim;
    struct reb_bridge* bridge;
};

struct roundtrip_result {
    const char* method_name;
    double time_residual_yr;
    double max_position_residual_au;
    double max_velocity_residual_au_per_yr;
};

static struct reb_particle particle_from_state(const struct galilean_body_state* state) {
    struct reb_particle particle;
    memset(&particle, 0, sizeof(particle));
    particle.m = state->mass_msun;
    particle.x = state->x_au;
    particle.y = state->y_au;
    particle.z = state->z_au;
    particle.vx = state->vx_au_yr;
    particle.vy = state->vy_au_yr;
    particle.vz = state->vz_au_yr;
    return particle;
}

static void add_cartesian(struct reb_simulation* sim, const struct reb_particle* particle) {
    reb_simulation_add_fmt(
        sim,
        "m x y z vx vy vz",
        particle->m,
        particle->x,
        particle->y,
        particle->z,
        particle->vx,
        particle->vy,
        particle->vz
    );
}

static struct reb_particle galilean_host_com(void) {
    struct reb_particle com;
    memset(&com, 0, sizeof(com));
    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        const struct reb_particle particle = particle_from_state(&GALILEAN_BODIES[i]);
        com.m += particle.m;
        com.x += particle.m * particle.x;
        com.y += particle.m * particle.y;
        com.z += particle.m * particle.z;
        com.vx += particle.m * particle.vx;
        com.vy += particle.m * particle.vy;
        com.vz += particle.m * particle.vz;
    }
    com.x /= com.m;
    com.y /= com.m;
    com.z /= com.m;
    com.vx /= com.m;
    com.vy /= com.m;
    com.vz /= com.m;
    return com;
}

static struct reb_particle subtract_particle(
    const struct reb_particle* particle,
    const struct reb_particle* origin
) {
    struct reb_particle result = *particle;
    result.x -= origin->x;
    result.y -= origin->y;
    result.z -= origin->z;
    result.vx -= origin->vx;
    result.vy -= origin->vy;
    result.vz -= origin->vz;
    return result;
}

static struct reb_particle add_particle(
    const struct reb_particle* particle,
    const struct reb_particle* origin
) {
    struct reb_particle result = *particle;
    result.x += origin->x;
    result.y += origin->y;
    result.z += origin->z;
    result.vx += origin->vx;
    result.vy += origin->vy;
    result.vz += origin->vz;
    return result;
}

static struct bridge_state make_bridge_state(enum reb_bridge_method method) {
    const double dt_outer = 0.5 / 365.25;
    const double dt_inner = 0.0125 / 365.25;
    struct bridge_state state;
    memset(&state, 0, sizeof(state));

    state.main_sim = reb_simulation_create();
    state.sub_sim = reb_simulation_create();
    if (!state.main_sim || !state.sub_sim) return state;
    state.main_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    state.sub_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(state.main_sim) != 0 ||
        reb_bridge_set_integrator_whfast(state.sub_sim) != 0) {
        return state;
    }
    state.main_sim->dt = dt_outer;
    state.sub_sim->dt = dt_inner;
    state.main_sim->exact_finish_time = 1;
    state.sub_sim->exact_finish_time = 1;

    struct reb_particle sun;
    memset(&sun, 0, sizeof(sun));
    sun.m = 1.0;
    const struct reb_particle host = galilean_host_com();
    add_cartesian(state.main_sim, &sun);
    add_cartesian(state.main_sim, &host);
    reb_simulation_move_to_com(state.main_sim);

    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        const struct reb_particle absolute = particle_from_state(&GALILEAN_BODIES[i]);
        const struct reb_particle relative = subtract_particle(&absolute, &host);
        add_cartesian(state.sub_sim, &relative);
    }
    reb_simulation_move_to_com(state.sub_sim);

    state.bridge = reb_bridge_create(state.main_sim, dt_outer);
    if (!state.bridge) return state;
    reb_bridge_set_owns_main(state.bridge, 1);
    if (reb_bridge_set_method(state.bridge, method) != 0 ||
        reb_bridge_add_subsystem(state.bridge, 1, state.sub_sim, dt_inner, 1) != 0) {
        return state;
    }
    return state;
}

static void free_bridge_state(struct bridge_state* state) {
    if (state->bridge) {
        reb_bridge_free(state->bridge);
    } else {
        if (state->main_sim) reb_simulation_free(state->main_sim);
        if (state->sub_sim) reb_simulation_free(state->sub_sim);
    }
    memset(state, 0, sizeof(*state));
}

static void fill_inertial_bodies(
    const struct bridge_state* state,
    struct reb_particle bodies[BODY_COUNT]
) {
    const struct reb_particle host = state->main_sim->particles[1];
    bodies[BODY_SUN] = state->main_sim->particles[0];
    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        bodies[i + 1] = add_particle(&state->sub_sim->particles[i], &host);
    }
}

static double position_delta(
    const struct reb_particle* initial,
    const struct reb_particle* final
) {
    const double dx = final->x - initial->x;
    const double dy = final->y - initial->y;
    const double dz = final->z - initial->z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

static double velocity_delta(
    const struct reb_particle* initial,
    const struct reb_particle* final
) {
    const double dvx = final->vx - initial->vx;
    const double dvy = final->vy - initial->vy;
    const double dvz = final->vz - initial->vz;
    return sqrt(dvx * dvx + dvy * dvy + dvz * dvz);
}

static int measure_roundtrip(
    enum reb_bridge_method method,
    const char* method_name,
    struct roundtrip_result* result
) {
    const int steps = 1000;
    const double dt = 0.5 / 365.25;
    struct bridge_state state = make_bridge_state(method);
    if (!state.bridge || !state.main_sim || !state.sub_sim) {
        free_bridge_state(&state);
        return -1;
    }

    reb_simulation_synchronize(state.main_sim);
    reb_simulation_synchronize(state.sub_sim);
    struct reb_particle initial[BODY_COUNT];
    fill_inertial_bodies(&state, initial);
    const double initial_main_time = state.main_sim->t;
    const double initial_sub_time = state.sub_sim->t;

    for (int i = 0; i < steps; i++) {
        if (reb_bridge_step(state.bridge, dt) != 0) {
            free_bridge_state(&state);
            return -1;
        }
    }
    for (int i = 0; i < steps; i++) {
        if (reb_bridge_step(state.bridge, -dt) != 0) {
            free_bridge_state(&state);
            return -1;
        }
    }

    reb_simulation_synchronize(state.main_sim);
    reb_simulation_synchronize(state.sub_sim);
    struct reb_particle final[BODY_COUNT];
    fill_inertial_bodies(&state, final);

    result->method_name = method_name;
    result->time_residual_yr = fmax(
        fabs(state.main_sim->t - initial_main_time),
        fabs(state.sub_sim->t - initial_sub_time)
    );
    result->max_position_residual_au = 0.0;
    result->max_velocity_residual_au_per_yr = 0.0;
    for (int i = 0; i < BODY_COUNT; i++) {
        const double dr = position_delta(&initial[i], &final[i]);
        const double dv = velocity_delta(&initial[i], &final[i]);
        if (dr > result->max_position_residual_au) result->max_position_residual_au = dr;
        if (dv > result->max_velocity_residual_au_per_yr) {
            result->max_velocity_residual_au_per_yr = dv;
        }
    }

    free_bridge_state(&state);
    if (!isfinite(result->time_residual_yr) ||
        !isfinite(result->max_position_residual_au) ||
        !isfinite(result->max_velocity_residual_au_per_yr)) {
        return -1;
    }
    if (result->time_residual_yr > 1.0e-9 ||
        result->max_position_residual_au > 1.0e-9 ||
        result->max_velocity_residual_au_per_yr > 1.0e-9) {
        return -1;
    }
    return 0;
}

static int parse_csv_path(int argc, char** argv, const char** csv_path) {
    if (argc == 3 && strcmp(argv[1], "--csv") == 0) {
        *csv_path = argv[2];
        return 0;
    }
    fprintf(stderr, "usage: %s --csv PATH\n", argv[0]);
    return -1;
}

int main(int argc, char** argv) {
    const char* csv_path = NULL;
    if (parse_csv_path(argc, argv, &csv_path) != 0) return EXIT_FAILURE;

    struct roundtrip_result results[2];
    if (measure_roundtrip(REB_BRIDGE_METHOD_ORDER2, "order2", &results[0]) != 0 ||
        measure_roundtrip(REB_BRIDGE_METHOD_EOS_LF4_2, "eos_lf4_2", &results[1]) != 0) {
        fprintf(stderr, "Galilean reversibility measurement failed\n");
        return EXIT_FAILURE;
    }

    FILE* csv = fopen(csv_path, "w");
    if (!csv) {
        fprintf(stderr, "failed to open CSV path: %s\n", csv_path);
        return EXIT_FAILURE;
    }
    fprintf(
        csv,
        "method,outer_steps_each_direction,dt_outer_days,dt_inner_days,"
        "time_residual_yr,max_position_residual_au,max_velocity_residual_au_per_yr\n"
    );
    for (int i = 0; i < 2; i++) {
        fprintf(
            csv,
            "%s,1000,0.5,0.0125,%.17g,%.17g,%.17g\n",
            results[i].method_name,
            results[i].time_residual_yr,
            results[i].max_position_residual_au,
            results[i].max_velocity_residual_au_per_yr
        );
        printf(
            "%s time=%.17g position=%.17g velocity=%.17g\n",
            results[i].method_name,
            results[i].time_residual_yr,
            results[i].max_position_residual_au,
            results[i].max_velocity_residual_au_per_yr
        );
    }
    fclose(csv);
    return EXIT_SUCCESS;
}
