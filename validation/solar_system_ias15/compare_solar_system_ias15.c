/* SPDX-License-Identifier: GPL-3.0-or-later */

/**
 * compare_solar_system_ias15.c - Bridge vs IAS15 validation for the solar-system demo.
 *
 * Produces a CSV with bridge and IAS15 diagnostics sampled at identical times.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "rebound_bridge.h"
#include "../common/val_cli.h"
#include "../common/val_particle.h"
#include "../common/solar_system_fixture.h"
#include "../validation_timer.h"


enum main_body {
    MAIN_SUN = 0,
    MAIN_MERCURY = 1,
    MAIN_VENUS = 2,
    MAIN_EARTH_MOON_HOST = 3,
    MAIN_MARS = 4,
    MAIN_JUPITER_HOST = 5,
    MAIN_SATURN = 6,
    MAIN_URANUS = 7,
    MAIN_NEPTUNE = 8
};

enum physical_body {
    BODY_SUN = 0,
    BODY_MERCURY = 1,
    BODY_VENUS = 2,
    BODY_EARTH = 3,
    BODY_MOON = 4,
    BODY_MARS = 5,
    BODY_JUPITER = 6,
    BODY_IO = 7,
    BODY_EUROPA = 8,
    BODY_GANYMEDE = 9,
    BODY_SATURN = 10,
    BODY_URANUS = 11,
    BODY_NEPTUNE = 12,
    BODY_COUNT = 13
};

enum earth_moon_body {
    EM_EARTH = 0,
    EM_MOON = 1
};

enum jovian_body {
    JOVIAN_JUPITER = 0,
    JOVIAN_IO = 1,
    JOVIAN_EUROPA = 2,
    JOVIAN_GANYMEDE = 3
};

struct bridge_state {
    struct reb_simulation* main_sim;
    struct reb_simulation* earth_moon;
    struct reb_simulation* jovian;
    struct reb_bridge* bridge;
};

struct body_metrics {
    double total_energy;
    double momentum_scale;
    struct reb_vec3d total_p;
    struct reb_vec3d total_l;
    struct reb_particle earth_moon_barycenter;
    struct reb_particle jovian_barycenter;
};

static double wrap_degrees(double angle) {
    angle = fmod(angle, 360.0);
    if (angle < 0.0) angle += 360.0;
    return angle;
}

static double centered_delta_degrees(double a, double b) {
    double delta = a - b;
    while (delta <= -180.0) delta += 360.0;
    while (delta > 180.0) delta -= 360.0;
    return delta;
}

static double distance_particles(const struct reb_particle* a, const struct reb_particle* b) {
    const double dx = a->x - b->x;
    const double dy = a->y - b->y;
    const double dz = a->z - b->z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

static double velocity_delta(const struct reb_particle* a, const struct reb_particle* b) {
    const double dvx = a->vx - b->vx;
    const double dvy = a->vy - b->vy;
    const double dvz = a->vz - b->vz;
    return sqrt(dvx * dvx + dvy * dvy + dvz * dvz);
}

static double relative_position_delta(
    const struct reb_particle* bridge_primary,
    const struct reb_particle* bridge_secondary,
    const struct reb_particle* ias15_primary,
    const struct reb_particle* ias15_secondary
) {
    const double dx =
        (bridge_secondary->x - bridge_primary->x) -
        (ias15_secondary->x - ias15_primary->x);
    const double dy =
        (bridge_secondary->y - bridge_primary->y) -
        (ias15_secondary->y - ias15_primary->y);
    const double dz =
        (bridge_secondary->z - bridge_primary->z) -
        (ias15_secondary->z - ias15_primary->z);
    return sqrt(dx * dx + dy * dy + dz * dz);
}


static void add_cartesian(
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


static struct reb_simulation* make_main_solar_system(void) {
    struct reb_simulation* sim = reb_simulation_create();

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    reb_simulation_add_fmt(sim, "m", 1.0);
    add_cartesian(sim, 1.660120e-7, -2.15201304377667713e-01, -4.09207615700272087e-01, -1.37032695919236600e-02, 7.02443404629697010e+00, -4.28716799788450942e+00, -9.94632701667513053e-01);
    add_cartesian(sim, 2.4478383e-6, 8.88772496163565545e-02, -7.21762382317080675e-01, -1.50440565750619305e-02, 7.28255910865529810e+00, 8.76733642322772511e-01, -4.08156573658943866e-01);
    add_cartesian(sim, 3.0034896e-6 + 3.6940e-8, -1.74258179548951841e-01, 9.67812101172537687e-01, -5.43185959007052695e-05, -6.28907787760335690e+00, -1.13619641904038438e+00, 1.15085854523419378e-04);
    add_cartesian(sim, 3.227151e-7, 3.40579676862215075e-01, -1.38700201594525407e+00, -3.74172267877010828e-02, 5.15711521620981550e+00, 1.65830921706158052e+00, -9.17089273335817273e-02);
    add_cartesian(sim, 9.547919e-4 + 2.388e-8 + 2.527e-8 + 7.804e-8, -1.69400240075433639e+00, 4.92888237470638568e+00, 1.74262370297428894e-02, -2.64080244023013400e+00, -7.67319806325283849e-01, 6.22885142647012779e-02);
    add_cartesian(sim, 2.8588567e-4, 9.50734336522187107e+00, 2.57738457802152210e-01, -3.82935528507753087e-01, -1.69141631745304527e-01, 2.03249650346815258e+00, -2.86682018818476998e-02);
    add_cartesian(sim, 4.366244e-5, 9.88031797519427535e+00, 1.68000154098179486e+01, -6.57197991373143220e-02, -1.25149298434788658e+00, 6.61240423623731699e-01, 1.86840216329132620e-02);
    add_cartesian(sim, 5.151389e-5, 2.98721200408889800e+01, 5.18939463241125942e-01, -6.99067175203629043e-01, -2.98764379233138963e-02, 1.15313346010349438e+00, -2.30704057616862752e-02);

    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_sub_earth_moon_horizons(double inner_mu_scale, unsigned int inner_whfast_corrector) {
    struct reb_simulation* sim = reb_simulation_create();

    sim->G = REB_BRIDGE_G_AU_YR_MSUN * inner_mu_scale;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->ri_whfast.corrector = inner_whfast_corrector;

    add_cartesian(sim, 3.0034896e-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    add_cartesian(sim, 3.6940e-8, 9.53035527126304600e-04, 2.17512105974572998e-03, 2.09676507944608496e-04, -2.09285374018137460e-01, 8.77025716574588027e-02, 1.16388188872971468e-03);

    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_sub_jovian_horizons(void) {
    struct reb_simulation* sim = reb_simulation_create();

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    reb_simulation_add_fmt(sim, "m", 9.547919e-4);
    add_cartesian(sim, 2.388e-8, 2.48120485765584205e-03, -1.34556438159649208e-03, -1.31630007347439604e-05, 1.75341752994026812e+00, 3.20163492363371871e+00, 1.39542430489110164e-01);
    add_cartesian(sim, 2.527e-8, 5.51813635549570501e-04, -4.47551321186325738e-03, -1.27610904592036212e-04, 2.85755306406471643e+00, 3.72439324002732552e-01, 7.16871253077594472e-02);
    add_cartesian(sim, 7.804e-8, 6.76397567656043287e-03, -2.29852901298470890e-03, 1.01873701583559195e-05, 7.40898140660848692e-01, 2.17447712709783314e+00, 9.37799108604124437e-02);

    reb_simulation_move_to_com(sim);
    return sim;
}

static struct bridge_state make_bridge_state(
    double dt_outer,
    double dt_earth_moon,
    double dt_jovian,
    double inner_mu_scale,
    unsigned int inner_whfast_corrector
) {
    struct bridge_state state;
    memset(&state, 0, sizeof(state));

    state.main_sim = make_main_solar_system();
    state.earth_moon = make_sub_earth_moon_horizons(inner_mu_scale, inner_whfast_corrector);
    state.jovian = make_sub_jovian_horizons();
    if (!state.main_sim || !state.earth_moon || !state.jovian) return state;

    state.bridge = reb_bridge_create(state.main_sim, dt_outer);
    if (!state.bridge) return state;

    reb_bridge_set_owns_main(state.bridge, 1);
    if (reb_bridge_add_subsystem(state.bridge, MAIN_EARTH_MOON_HOST, state.earth_moon, dt_earth_moon, 1) != 0) return state;
    if (reb_bridge_add_subsystem(state.bridge, MAIN_JUPITER_HOST, state.jovian, dt_jovian, 1) != 0) return state;
    return state;
}

static void free_bridge_state(struct bridge_state* state) {
    if (!state) return;
    if (state->bridge) {
        reb_bridge_free(state->bridge);
    } else {
        if (state->main_sim) reb_simulation_free(state->main_sim);
        if (state->earth_moon) reb_simulation_free(state->earth_moon);
        if (state->jovian) reb_simulation_free(state->jovian);
    }
    memset(state, 0, sizeof(*state));
}

static struct reb_simulation* make_ias15_reference(void) {
    struct reb_simulation* sim = reb_simulation_create();

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->exact_finish_time = 1;

    struct reb_particle bodies[BODY_COUNT];
    ss_populate_bodies(bodies);
    for (int i = 0; i < BODY_COUNT; i++) {
        reb_simulation_add(sim, bodies[i]);
    }

    reb_simulation_move_to_com(sim);
    return sim;
}


static void fill_bridge_bodies(const struct bridge_state* state, struct reb_particle bodies[BODY_COUNT]) {
    const struct reb_simulation* main_sim = state->main_sim;
    const struct reb_simulation* earth_moon = state->earth_moon;
    const struct reb_simulation* jovian = state->jovian;
    const struct reb_particle* host_em = &main_sim->particles[MAIN_EARTH_MOON_HOST];
    const struct reb_particle* host_jov = &main_sim->particles[MAIN_JUPITER_HOST];

    bodies[BODY_SUN] = main_sim->particles[MAIN_SUN];
    bodies[BODY_MERCURY] = main_sim->particles[MAIN_MERCURY];
    bodies[BODY_VENUS] = main_sim->particles[MAIN_VENUS];
    bodies[BODY_EARTH] = val_particle_add(host_em, &earth_moon->particles[EM_EARTH]);
    bodies[BODY_MOON] = val_particle_add(host_em, &earth_moon->particles[EM_MOON]);
    bodies[BODY_MARS] = main_sim->particles[MAIN_MARS];
    bodies[BODY_JUPITER] = val_particle_add(host_jov, &jovian->particles[JOVIAN_JUPITER]);
    bodies[BODY_IO] = val_particle_add(host_jov, &jovian->particles[JOVIAN_IO]);
    bodies[BODY_EUROPA] = val_particle_add(host_jov, &jovian->particles[JOVIAN_EUROPA]);
    bodies[BODY_GANYMEDE] = val_particle_add(host_jov, &jovian->particles[JOVIAN_GANYMEDE]);
    bodies[BODY_SATURN] = main_sim->particles[MAIN_SATURN];
    bodies[BODY_URANUS] = main_sim->particles[MAIN_URANUS];
    bodies[BODY_NEPTUNE] = main_sim->particles[MAIN_NEPTUNE];
}

static void fill_reference_bodies(const struct reb_simulation* sim, struct reb_particle bodies[BODY_COUNT]) {
    for (int i = 0; i < BODY_COUNT; i++) {
        bodies[i] = sim->particles[i];
    }
}

static struct body_metrics compute_metrics(const struct reb_particle bodies[BODY_COUNT], double G) {
    struct body_metrics metrics;
    memset(&metrics, 0, sizeof(metrics));

    for (int i = 0; i < BODY_COUNT; i++) {
        const struct reb_particle* p = &bodies[i];
        const double v2 = p->vx * p->vx + p->vy * p->vy + p->vz * p->vz;
        metrics.total_energy += 0.5 * p->m * v2;
        metrics.total_p.x += p->m * p->vx;
        metrics.total_p.y += p->m * p->vy;
        metrics.total_p.z += p->m * p->vz;
        metrics.momentum_scale += p->m * sqrt(v2);
        metrics.total_l.x += p->m * (p->y * p->vz - p->z * p->vy);
        metrics.total_l.y += p->m * (p->z * p->vx - p->x * p->vz);
        metrics.total_l.z += p->m * (p->x * p->vy - p->y * p->vx);
    }

    for (int i = 0; i < BODY_COUNT; i++) {
        for (int j = i + 1; j < BODY_COUNT; j++) {
            metrics.total_energy -= G * bodies[i].m * bodies[j].m / distance_particles(&bodies[i], &bodies[j]);
        }
    }

    const int em_ids[] = { BODY_EARTH, BODY_MOON };
    const int jov_ids[] = { BODY_JUPITER, BODY_IO, BODY_EUROPA, BODY_GANYMEDE };
    metrics.earth_moon_barycenter = val_barycenter(bodies, em_ids, 2);
    metrics.jovian_barycenter = val_barycenter(bodies, jov_ids, 4);
    return metrics;
}

static double relative_energy_error(double energy, double energy0) {
    const double denom = fabs(energy0) > 0.0 ? fabs(energy0) : 1.0;
    return (energy - energy0) / denom;
}

static double relative_l_error(struct reb_vec3d l, struct reb_vec3d l0) {
    const double denom = sqrt(l0.x * l0.x + l0.y * l0.y + l0.z * l0.z);
    if (denom <= 0.0) return 0.0;
    const double dx = l.x - l0.x;
    const double dy = l.y - l0.y;
    const double dz = l.z - l0.z;
    return sqrt(dx * dx + dy * dy + dz * dz) / denom;
}

static double relative_momentum_error(
    struct reb_vec3d p,
    struct reb_vec3d p0,
    double momentum_scale0
) {
    const double denom = momentum_scale0 > 0.0 ? momentum_scale0 : 1.0;
    const double dx = p.x - p0.x;
    const double dy = p.y - p0.y;
    const double dz = p.z - p0.z;
    return sqrt(dx * dx + dy * dy + dz * dz) / denom;
}

static double laplace_angle_from_bodies(const struct reb_particle bodies[BODY_COUNT]) {
    const struct reb_orbit io = reb_orbit_from_particle(REB_BRIDGE_G_AU_YR_MSUN, bodies[BODY_IO], bodies[BODY_JUPITER]);
    const struct reb_orbit europa = reb_orbit_from_particle(REB_BRIDGE_G_AU_YR_MSUN, bodies[BODY_EUROPA], bodies[BODY_JUPITER]);
    const struct reb_orbit ganymede = reb_orbit_from_particle(REB_BRIDGE_G_AU_YR_MSUN, bodies[BODY_GANYMEDE], bodies[BODY_JUPITER]);
    return wrap_degrees((io.l - 3.0 * europa.l + 2.0 * ganymede.l) * 180.0 / M_PI);
}

static double heliocentric_radius(const struct reb_particle* body, const struct reb_particle* sun) {
    return distance_particles(body, sun);
}

static double heliocentric_phase_degrees(const struct reb_particle* body, const struct reb_particle* sun) {
    return wrap_degrees(atan2(body->y - sun->y, body->x - sun->x) * 180.0 / M_PI);
}

static const char* physical_body_name(int body_index) {
    switch (body_index) {
        case BODY_SUN:
            return "Sun";
        case BODY_MERCURY:
            return "Mercury";
        case BODY_VENUS:
            return "Venus";
        case BODY_EARTH:
            return "Earth";
        case BODY_MOON:
            return "Moon";
        case BODY_MARS:
            return "Mars";
        case BODY_JUPITER:
            return "Jupiter";
        case BODY_IO:
            return "Io";
        case BODY_EUROPA:
            return "Europa";
        case BODY_GANYMEDE:
            return "Ganymede";
        case BODY_SATURN:
            return "Saturn";
        case BODY_URANUS:
            return "Uranus";
        case BODY_NEPTUNE:
            return "Neptune";
    }
    return "Unknown";
}

static int write_initial_state_csv(
    const char* out_path,
    const struct reb_particle bridge_bodies[BODY_COUNT],
    const struct reb_particle ias15_bodies[BODY_COUNT],
    double* max_abs_mass_diff,
    double* max_dr_au,
    double* max_dv_au_per_yr
) {
    FILE* fp = fopen(out_path, "w");
    if (!fp) {
        fprintf(stderr, "failed to open initial-state CSV: %s\n", out_path);
        return -1;
    }

    *max_abs_mass_diff = 0.0;
    *max_dr_au = 0.0;
    *max_dv_au_per_yr = 0.0;
    fprintf(
        fp,
        "body,bridge_m,ias15_m,mass_diff,"
        "bridge_x,bridge_y,bridge_z,bridge_vx,bridge_vy,bridge_vz,"
        "ias15_x,ias15_y,ias15_z,ias15_vx,ias15_vy,ias15_vz,"
        "dr_au,dv_au_per_yr\n"
    );

    for (int i = 0; i < BODY_COUNT; i++) {
        const struct reb_particle* bridge = &bridge_bodies[i];
        const struct reb_particle* ias15 = &ias15_bodies[i];
        const double mass_diff = bridge->m - ias15->m;
        const double dr = distance_particles(bridge, ias15);
        const double dv = velocity_delta(bridge, ias15);
        if (fabs(mass_diff) > *max_abs_mass_diff) *max_abs_mass_diff = fabs(mass_diff);
        if (dr > *max_dr_au) *max_dr_au = dr;
        if (dv > *max_dv_au_per_yr) *max_dv_au_per_yr = dv;
        fprintf(
            fp,
            "%s,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g\n",
            physical_body_name(i),
            bridge->m,
            ias15->m,
            mass_diff,
            bridge->x,
            bridge->y,
            bridge->z,
            bridge->vx,
            bridge->vy,
            bridge->vz,
            ias15->x,
            ias15->y,
            ias15->z,
            ias15->vx,
            ias15->vy,
            ias15->vz,
            dr,
            dv
        );
    }

    fclose(fp);
    return 0;
}

static void write_csv_header(FILE* fp) {
    fprintf(fp,
        "bridge_cross_kick,"
        "force_gradient_scale,"
        "bridge_corrector,"
        "bridge_corrector_scale,"
        "inner_mu_scale,"
        "inner_whfast_corrector,"
        "time_yr,"
        "bridge_step_seconds,ias15_step_seconds,step_speedup,"
        "bridge_cumulative_seconds,ias15_cumulative_seconds,cumulative_speedup,"
        "bridge_energy_rel,ias15_energy_rel,energy_rel_diff,"
        "bridge_l_rel,ias15_l_rel,l_rel_diff,"
        "bridge_p_rel,ias15_p_rel,p_rel_diff,"
        "earth_abs_error_au,moon_abs_error_au,moon_relative_error_au,"
        "bridge_em_distance_au,ias15_em_distance_au,em_distance_diff_au,"
        "bridge_em_a_au,ias15_em_a_au,em_a_diff_au,"
        "bridge_em_e,ias15_em_e,em_e_diff,"
        "bridge_laplace_deg,ias15_laplace_deg,laplace_diff_deg,"
        "earth_barycenter_error_au,jupiter_barycenter_error_au,"
        "earth_phase_diff_deg,jupiter_phase_diff_deg,"
        "earth_radius_diff_au,jupiter_radius_diff_au\n"
    );
}

static void write_csv_row(
    FILE* fp,
    const char* bridge_cross_kick,
    double force_gradient_scale,
    const char* bridge_corrector,
    double bridge_corrector_scale,
    double inner_mu_scale,
    unsigned int inner_whfast_corrector,
    double time_yr,
    double bridge_step_seconds,
    double ias15_step_seconds,
    double bridge_cumulative_seconds,
    double ias15_cumulative_seconds,
    double bridge_energy_rel,
    double ias15_energy_rel,
    double bridge_l_rel,
    double ias15_l_rel,
    double bridge_p_rel,
    double ias15_p_rel,
    double earth_abs_error_au,
    double moon_abs_error_au,
    double moon_relative_error_au,
    double bridge_em_distance,
    double ias15_em_distance,
    double bridge_em_a_au,
    double ias15_em_a_au,
    double bridge_em_e,
    double ias15_em_e,
    double bridge_laplace_deg,
    double ias15_laplace_deg,
    double earth_barycenter_error_au,
    double jupiter_barycenter_error_au,
    double earth_phase_diff_deg,
    double jupiter_phase_diff_deg,
    double earth_radius_diff_au,
    double jupiter_radius_diff_au
) {
    const double step_speedup = bridge_step_seconds > 0.0 ? ias15_step_seconds / bridge_step_seconds : 0.0;
    const double cumulative_speedup = bridge_cumulative_seconds > 0.0 ? ias15_cumulative_seconds / bridge_cumulative_seconds : 0.0;

#define WRITE_CSV_VALUE(value) fprintf(fp, ",%.16e", (value))
    fprintf(fp, "%s", bridge_cross_kick);
    WRITE_CSV_VALUE(force_gradient_scale);
    fprintf(fp, ",%s", bridge_corrector);
    WRITE_CSV_VALUE(bridge_corrector_scale);
    WRITE_CSV_VALUE(inner_mu_scale);
    fprintf(fp, ",%u", inner_whfast_corrector);
    WRITE_CSV_VALUE(time_yr);
    WRITE_CSV_VALUE(bridge_step_seconds);
    WRITE_CSV_VALUE(ias15_step_seconds);
    WRITE_CSV_VALUE(step_speedup);
    WRITE_CSV_VALUE(bridge_cumulative_seconds);
    WRITE_CSV_VALUE(ias15_cumulative_seconds);
    WRITE_CSV_VALUE(cumulative_speedup);
    WRITE_CSV_VALUE(bridge_energy_rel);
    WRITE_CSV_VALUE(ias15_energy_rel);
    WRITE_CSV_VALUE(bridge_energy_rel - ias15_energy_rel);
    WRITE_CSV_VALUE(bridge_l_rel);
    WRITE_CSV_VALUE(ias15_l_rel);
    WRITE_CSV_VALUE(bridge_l_rel - ias15_l_rel);
    WRITE_CSV_VALUE(bridge_p_rel);
    WRITE_CSV_VALUE(ias15_p_rel);
    WRITE_CSV_VALUE(bridge_p_rel - ias15_p_rel);
    WRITE_CSV_VALUE(earth_abs_error_au);
    WRITE_CSV_VALUE(moon_abs_error_au);
    WRITE_CSV_VALUE(moon_relative_error_au);
    WRITE_CSV_VALUE(bridge_em_distance);
    WRITE_CSV_VALUE(ias15_em_distance);
    WRITE_CSV_VALUE(bridge_em_distance - ias15_em_distance);
    WRITE_CSV_VALUE(bridge_em_a_au);
    WRITE_CSV_VALUE(ias15_em_a_au);
    WRITE_CSV_VALUE(bridge_em_a_au - ias15_em_a_au);
    WRITE_CSV_VALUE(bridge_em_e);
    WRITE_CSV_VALUE(ias15_em_e);
    WRITE_CSV_VALUE(bridge_em_e - ias15_em_e);
    WRITE_CSV_VALUE(bridge_laplace_deg);
    WRITE_CSV_VALUE(ias15_laplace_deg);
    WRITE_CSV_VALUE(bridge_laplace_deg - ias15_laplace_deg);
    WRITE_CSV_VALUE(earth_barycenter_error_au);
    WRITE_CSV_VALUE(jupiter_barycenter_error_au);
    WRITE_CSV_VALUE(earth_phase_diff_deg);
    WRITE_CSV_VALUE(jupiter_phase_diff_deg);
    WRITE_CSV_VALUE(earth_radius_diff_au);
    WRITE_CSV_VALUE(jupiter_radius_diff_au);
    fprintf(fp, "\n");
#undef WRITE_CSV_VALUE
}

static void print_usage(const char* argv0) {
    fprintf(
        stderr,
        "usage: %s [--years N] [--samples N] [--dt-outer-days N] [--earth-moon-ratio N] [--jovian-ratio N] [--bridge-method order2|eos-lf4-2|eos-lf8-6-4|forward-fg4] [--bridge-cross-kick plain|force-gradient] [--force-gradient-scale N] [--bridge-corrector none|commutator-ab] [--bridge-corrector-scale N] [--inner-mu-scale N] [--inner-whfast-corrector 0|11] [--csv PATH] [--initial-state-csv PATH] [--force-gradient-check] [--force-gradient-check-epsilon N]\n",
        argv0
    );
}


static const char* cross_kick_kernel_name(enum reb_bridge_cross_kick_kernel kernel) {
    switch (kernel) {
        case REB_BRIDGE_CROSS_KICK_PLAIN:
            return "plain";
        case REB_BRIDGE_CROSS_KICK_FORCE_GRADIENT:
            return "force-gradient";
        case REB_BRIDGE_CROSS_KICK_INVALID:
            return "invalid";
    }
    return "unknown";
}

static const char* bridge_corrector_name(enum reb_bridge_corrector corrector) {
    switch (corrector) {
        case REB_BRIDGE_CORRECTOR_NONE:
            return "none";
        case REB_BRIDGE_CORRECTOR_COMMUTATOR_AB:
            return "commutator-ab";
        case REB_BRIDGE_CORRECTOR_INVALID:
            return "invalid";
    }
    return "unknown";
}

int main(int argc, char** argv) {
    double years = 2000.0;
    int samples = 2000;
    double dt_outer_days = 1.0;
    const char* out_path = "validation/solar_system_ias15/out/solar_system_ias15_compare.csv";
    int earth_moon_ratio = 20;
    int jovian_ratio = 50;
    double inner_mu_scale = 1.0;
    unsigned int inner_whfast_corrector = 0;
    const char* initial_state_path = NULL;
    int force_gradient_check = 0;
    double force_gradient_check_epsilon = 1.0e-7;
    enum reb_bridge_method method = REB_BRIDGE_METHOD_ORDER2;
    enum reb_bridge_cross_kick_kernel cross_kick_kernel = REB_BRIDGE_CROSS_KICK_PLAIN;
    double force_gradient_scale = 1.0;
    enum reb_bridge_corrector bridge_corrector = REB_BRIDGE_CORRECTOR_NONE;
    double bridge_corrector_scale = 1.0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--years") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            years = atof(argv[++i]);
            if (years <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--samples") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            samples = atoi(argv[++i]);
            if (samples <= 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--dt-outer-days") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            dt_outer_days = atof(argv[++i]);
            if (dt_outer_days <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--earth-moon-ratio") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            earth_moon_ratio = atoi(argv[++i]);
            if (earth_moon_ratio <= 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--jovian-ratio") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            jovian_ratio = atoi(argv[++i]);
            if (jovian_ratio <= 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--bridge-method") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            const char* value = argv[++i];
            if (strcmp(value, "order2") == 0) {
                method = REB_BRIDGE_METHOD_ORDER2;
            } else if (strcmp(value, "eos-lf4-2") == 0) {
                method = REB_BRIDGE_METHOD_EOS_LF4_2;
            } else if (strcmp(value, "eos-lf8-6-4") == 0) {
                method = REB_BRIDGE_METHOD_EOS_LF8_6_4;
            } else if (strcmp(value, "forward-fg4") == 0) {
                method = REB_BRIDGE_METHOD_FORWARD_FG4;
            } else {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--bridge-cross-kick") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            const char* value = argv[++i];
            if (strcmp(value, "plain") == 0) {
                cross_kick_kernel = REB_BRIDGE_CROSS_KICK_PLAIN;
            } else if (strcmp(value, "force-gradient") == 0) {
                cross_kick_kernel = REB_BRIDGE_CROSS_KICK_FORCE_GRADIENT;
            } else {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--force-gradient-scale") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            force_gradient_scale = atof(argv[++i]);
            if (!isfinite(force_gradient_scale)) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--bridge-corrector") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            const char* value = argv[++i];
            if (strcmp(value, "none") == 0) {
                bridge_corrector = REB_BRIDGE_CORRECTOR_NONE;
            } else if (strcmp(value, "commutator-ab") == 0) {
                bridge_corrector = REB_BRIDGE_CORRECTOR_COMMUTATOR_AB;
            } else {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--bridge-corrector-scale") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            bridge_corrector_scale = atof(argv[++i]);
            if (!isfinite(bridge_corrector_scale)) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--inner-mu-scale") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            inner_mu_scale = atof(argv[++i]);
            if (inner_mu_scale <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--inner-whfast-corrector") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            const int corrector_value = atoi(argv[++i]);
            if (corrector_value != 0 && corrector_value != 11) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            inner_whfast_corrector = (unsigned int)corrector_value;
        } else if (strcmp(argv[i], "--csv") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            out_path = argv[++i];
        } else if (strcmp(argv[i], "--initial-state-csv") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            initial_state_path = argv[++i];
        } else if (strcmp(argv[i], "--force-gradient-check") == 0) {
            force_gradient_check = 1;
        } else if (strcmp(argv[i], "--force-gradient-check-epsilon") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            force_gradient_check_epsilon = atof(argv[++i]);
            if (force_gradient_check_epsilon <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    const double dt_outer = dt_outer_days / 365.25;
    const double dt_earth_moon = dt_outer / (double)earth_moon_ratio;
    const double dt_jovian = dt_outer / (double)jovian_ratio;

    struct bridge_state bridge_state = make_bridge_state(
        dt_outer,
        dt_earth_moon,
        dt_jovian,
        inner_mu_scale,
        inner_whfast_corrector
    );
    if (!bridge_state.bridge || !bridge_state.main_sim || !bridge_state.earth_moon || !bridge_state.jovian) {
        fprintf(stderr, "failed to create bridge state\n");
        free_bridge_state(&bridge_state);
        return EXIT_FAILURE;
    }
    if (reb_bridge_set_method(bridge_state.bridge, method) != 0) {
        fprintf(stderr, "failed to configure bridge method\n");
        free_bridge_state(&bridge_state);
        return EXIT_FAILURE;
    }
    if (reb_bridge_set_cross_kick_kernel(bridge_state.bridge, cross_kick_kernel) != 0) {
        fprintf(stderr, "failed to configure bridge cross-kick kernel\n");
        free_bridge_state(&bridge_state);
        return EXIT_FAILURE;
    }
    if (reb_bridge_set_force_gradient_scale(bridge_state.bridge, force_gradient_scale) != 0) {
        fprintf(stderr, "failed to configure bridge force-gradient scale\n");
        free_bridge_state(&bridge_state);
        return EXIT_FAILURE;
    }
    if (reb_bridge_set_corrector(bridge_state.bridge, bridge_corrector) != 0 ||
        reb_bridge_set_corrector_scale(bridge_state.bridge, bridge_corrector_scale) != 0) {
        fprintf(stderr, "failed to configure bridge corrector\n");
        free_bridge_state(&bridge_state);
        return EXIT_FAILURE;
    }

    struct reb_simulation* ias15 = make_ias15_reference();
    if (!ias15) {
        fprintf(stderr, "failed to create IAS15 reference state\n");
        free_bridge_state(&bridge_state);
        return EXIT_FAILURE;
    }

    if (initial_state_path) {
        struct reb_particle bridge_bodies[BODY_COUNT];
        struct reb_particle ias15_bodies[BODY_COUNT];
        fill_bridge_bodies(&bridge_state, bridge_bodies);
        fill_reference_bodies(ias15, ias15_bodies);

        double max_abs_mass_diff = 0.0;
        double max_dr_au = 0.0;
        double max_dv_au_per_yr = 0.0;
        if (write_initial_state_csv(
                initial_state_path,
                bridge_bodies,
                ias15_bodies,
                &max_abs_mass_diff,
                &max_dr_au,
                &max_dv_au_per_yr
            ) != 0) {
            reb_simulation_free(ias15);
            free_bridge_state(&bridge_state);
            return EXIT_FAILURE;
        }
        printf(
            "initial_state_csv=%s\nmax_abs_mass_diff=%.17g\nmax_dr_au=%.17g\nmax_dv_au_per_yr=%.17g\n",
            initial_state_path,
            max_abs_mass_diff,
            max_dr_au,
            max_dv_au_per_yr
        );
        reb_simulation_free(ias15);
        free_bridge_state(&bridge_state);
        return EXIT_SUCCESS;
    }

    if (force_gradient_check) {
        struct reb_bridge_force_gradient_check check;
        if (reb_bridge_force_gradient_finite_difference_check(
                bridge_state.bridge,
                0,
                force_gradient_check_epsilon,
                &check
            ) != 0) {
            reb_simulation_free(ias15);
            free_bridge_state(&bridge_state);
            return EXIT_FAILURE;
        }
        printf("force_gradient_check_epsilon=%.17g\n", force_gradient_check_epsilon);
        printf("force_gradient_max_abs_difference=%.17g\n", check.max_abs_difference);
        printf("force_gradient_max_relative_difference=%.17g\n", check.max_relative_difference);
        printf("force_gradient_max_abs_analytic=%.17g\n", check.max_abs_analytic);
        printf("force_gradient_max_abs_finite_difference=%.17g\n", check.max_abs_finite_difference);
        printf("force_gradient_internal_weighted_acc_residual=%.17g\n", check.force_gradient_internal_weighted_acc_residual);
        printf("force_gradient_main_weighted_acc_residual=%.17g\n", check.force_gradient_main_weighted_acc_residual);
        printf("cross_kick_internal_weighted_acc_residual=%.17g\n", check.internal_weighted_acc_residual);
        printf("cross_kick_main_weighted_acc_residual=%.17g\n", check.main_weighted_acc_residual);
        reb_simulation_free(ias15);
        free_bridge_state(&bridge_state);
        return EXIT_SUCCESS;
    }

    FILE* fp = fopen(out_path, "w");
    if (!fp) {
        fprintf(stderr, "failed to open output CSV: %s\n", out_path);
        reb_simulation_free(ias15);
        free_bridge_state(&bridge_state);
        return EXIT_FAILURE;
    }

    write_csv_header(fp);

    struct reb_particle bridge_initial[BODY_COUNT];
    struct reb_particle ias15_initial[BODY_COUNT];
    fill_bridge_bodies(&bridge_state, bridge_initial);
    fill_reference_bodies(ias15, ias15_initial);

    const struct body_metrics bridge0 = compute_metrics(bridge_initial, REB_BRIDGE_G_AU_YR_MSUN);
    const struct body_metrics ias150 = compute_metrics(ias15_initial, REB_BRIDGE_G_AU_YR_MSUN);
    const double bridge_energy0 = bridge0.total_energy;
    const double ias15_energy0 = ias150.total_energy;
    const struct reb_vec3d bridge_p0 = bridge0.total_p;
    const struct reb_vec3d ias15_p0 = ias150.total_p;
    const double bridge_p_scale0 = bridge0.momentum_scale;
    const double ias15_p_scale0 = ias150.momentum_scale;
    const struct reb_vec3d bridge_l0 = bridge0.total_l;
    const struct reb_vec3d ias15_l0 = ias150.total_l;
    double bridge_cumulative_seconds = 0.0;
    double ias15_cumulative_seconds = 0.0;

    if (bridge_corrector != REB_BRIDGE_CORRECTOR_NONE &&
        reb_bridge_apply_corrector(bridge_state.bridge, REB_BRIDGE_CORRECTOR_FORWARD) != 0) {
        fprintf(stderr, "bridge corrector warm-start failed\n");
        fclose(fp);
        reb_simulation_free(ias15);
        free_bridge_state(&bridge_state);
        return EXIT_FAILURE;
    }

    for (int i = 0; i < samples; i++) {
        const double target = years * (double)(i + 1) / (double)samples;
        const double bridge_start = monotonic_seconds();

        const double bridge_duration = target - bridge_state.main_sim->t;
        if (reb_bridge_advance(bridge_state.bridge, bridge_duration) != 0) {
            fprintf(stderr, "bridge advance failed at sample %d\n", i + 1);
            fclose(fp);
            reb_simulation_free(ias15);
            free_bridge_state(&bridge_state);
            return EXIT_FAILURE;
        }
        const double bridge_step_seconds = monotonic_seconds() - bridge_start;
        bridge_cumulative_seconds += bridge_step_seconds;
        const double ias15_start = monotonic_seconds();
        if (reb_simulation_integrate(ias15, target) != REB_STATUS_SUCCESS) {
            fprintf(stderr, "IAS15 integrate failed at sample %d\n", i + 1);
            fclose(fp);
            reb_simulation_free(ias15);
            free_bridge_state(&bridge_state);
            return EXIT_FAILURE;
        }
        const double ias15_step_seconds = monotonic_seconds() - ias15_start;
        ias15_cumulative_seconds += ias15_step_seconds;

        if (bridge_corrector != REB_BRIDGE_CORRECTOR_NONE &&
            reb_bridge_apply_corrector(bridge_state.bridge, REB_BRIDGE_CORRECTOR_INVERSE) != 0) {
            fprintf(stderr, "bridge corrector cold-output failed at sample %d\n", i + 1);
            fclose(fp);
            reb_simulation_free(ias15);
            free_bridge_state(&bridge_state);
            return EXIT_FAILURE;
        }

        struct reb_particle bridge_bodies[BODY_COUNT];
        struct reb_particle ias15_bodies[BODY_COUNT];
        fill_bridge_bodies(&bridge_state, bridge_bodies);
        fill_reference_bodies(ias15, ias15_bodies);

        const struct body_metrics bridge_metrics = compute_metrics(bridge_bodies, REB_BRIDGE_G_AU_YR_MSUN);
        const struct body_metrics ias15_metrics = compute_metrics(ias15_bodies, REB_BRIDGE_G_AU_YR_MSUN);

        const double bridge_energy_rel = relative_energy_error(bridge_metrics.total_energy, bridge_energy0);
        const double ias15_energy_rel = relative_energy_error(ias15_metrics.total_energy, ias15_energy0);
        const double bridge_l_rel = relative_l_error(bridge_metrics.total_l, bridge_l0);
        const double ias15_l_rel = relative_l_error(ias15_metrics.total_l, ias15_l0);
        const double bridge_p_rel = relative_momentum_error(bridge_metrics.total_p, bridge_p0, bridge_p_scale0);
        const double ias15_p_rel = relative_momentum_error(ias15_metrics.total_p, ias15_p0, ias15_p_scale0);
        const double earth_abs_error_au = distance_particles(&bridge_bodies[BODY_EARTH], &ias15_bodies[BODY_EARTH]);
        const double moon_abs_error_au = distance_particles(&bridge_bodies[BODY_MOON], &ias15_bodies[BODY_MOON]);
        const double moon_relative_error_au = relative_position_delta(
            &bridge_bodies[BODY_EARTH],
            &bridge_bodies[BODY_MOON],
            &ias15_bodies[BODY_EARTH],
            &ias15_bodies[BODY_MOON]
        );
        const double bridge_em_distance = distance_particles(&bridge_bodies[BODY_EARTH], &bridge_bodies[BODY_MOON]);
        const double ias15_em_distance = distance_particles(&ias15_bodies[BODY_EARTH], &ias15_bodies[BODY_MOON]);
        const struct reb_orbit bridge_em_orbit = reb_orbit_from_particle(
            REB_BRIDGE_G_AU_YR_MSUN,
            bridge_bodies[BODY_MOON],
            bridge_bodies[BODY_EARTH]
        );
        const struct reb_orbit ias15_em_orbit = reb_orbit_from_particle(
            REB_BRIDGE_G_AU_YR_MSUN,
            ias15_bodies[BODY_MOON],
            ias15_bodies[BODY_EARTH]
        );
        const double bridge_laplace_deg = laplace_angle_from_bodies(bridge_bodies);
        const double ias15_laplace_deg = laplace_angle_from_bodies(ias15_bodies);
        const double earth_barycenter_error_au = distance_particles(
            &bridge_state.main_sim->particles[MAIN_EARTH_MOON_HOST],
            &ias15_metrics.earth_moon_barycenter
        );
        const double jupiter_barycenter_error_au = distance_particles(
            &bridge_state.main_sim->particles[MAIN_JUPITER_HOST],
            &ias15_metrics.jovian_barycenter
        );
        const double earth_phase_diff_deg = centered_delta_degrees(
            heliocentric_phase_degrees(&bridge_metrics.earth_moon_barycenter, &bridge_bodies[BODY_SUN]),
            heliocentric_phase_degrees(&ias15_metrics.earth_moon_barycenter, &ias15_bodies[BODY_SUN])
        );
        const double jupiter_phase_diff_deg = centered_delta_degrees(
            heliocentric_phase_degrees(&bridge_metrics.jovian_barycenter, &bridge_bodies[BODY_SUN]),
            heliocentric_phase_degrees(&ias15_metrics.jovian_barycenter, &ias15_bodies[BODY_SUN])
        );
        const double earth_radius_diff_au =
            heliocentric_radius(&bridge_metrics.earth_moon_barycenter, &bridge_bodies[BODY_SUN]) -
            heliocentric_radius(&ias15_metrics.earth_moon_barycenter, &ias15_bodies[BODY_SUN]);
        const double jupiter_radius_diff_au =
            heliocentric_radius(&bridge_metrics.jovian_barycenter, &bridge_bodies[BODY_SUN]) -
            heliocentric_radius(&ias15_metrics.jovian_barycenter, &ias15_bodies[BODY_SUN]);

        write_csv_row(
            fp,
            cross_kick_kernel_name(cross_kick_kernel),
            force_gradient_scale,
            bridge_corrector_name(bridge_corrector),
            bridge_corrector_scale,
            inner_mu_scale,
            inner_whfast_corrector,
            target,
            bridge_step_seconds,
            ias15_step_seconds,
            bridge_cumulative_seconds,
            ias15_cumulative_seconds,
            bridge_energy_rel,
            ias15_energy_rel,
            bridge_l_rel,
            ias15_l_rel,
            bridge_p_rel,
            ias15_p_rel,
            earth_abs_error_au,
            moon_abs_error_au,
            moon_relative_error_au,
            bridge_em_distance,
            ias15_em_distance,
            bridge_em_orbit.a,
            ias15_em_orbit.a,
            bridge_em_orbit.e,
            ias15_em_orbit.e,
            bridge_laplace_deg,
            ias15_laplace_deg,
            earth_barycenter_error_au,
            jupiter_barycenter_error_au,
            earth_phase_diff_deg,
            jupiter_phase_diff_deg,
            earth_radius_diff_au,
            jupiter_radius_diff_au
        );

        if (bridge_corrector != REB_BRIDGE_CORRECTOR_NONE &&
            reb_bridge_apply_corrector(bridge_state.bridge, REB_BRIDGE_CORRECTOR_FORWARD) != 0) {
            fprintf(stderr, "bridge corrector re-warm failed at sample %d\n", i + 1);
            fclose(fp);
            reb_simulation_free(ias15);
            free_bridge_state(&bridge_state);
            return EXIT_FAILURE;
        }
    }

    fclose(fp);
    printf(
        "dt_outer_days=%.9f\nearth_moon_ratio=%d\njovian_ratio=%d\nbridge_method=%s\nbridge_cross_kick=%s\nforce_gradient_scale=%.17g\nbridge_corrector=%s\nbridge_corrector_scale=%.17g\ninner_mu_scale=%.17g\ninner_whfast_corrector=%u\nbridge_total_seconds=%.9f\nias15_total_seconds=%.9f\noverall_speedup=%.9f\n",
        dt_outer_days,
        earth_moon_ratio,
        jovian_ratio,
        val_bridge_method_name(method),
        cross_kick_kernel_name(cross_kick_kernel),
        force_gradient_scale,
        bridge_corrector_name(bridge_corrector),
        bridge_corrector_scale,
        inner_mu_scale,
        inner_whfast_corrector,
        bridge_cumulative_seconds,
        ias15_cumulative_seconds,
        bridge_cumulative_seconds > 0.0 ? ias15_cumulative_seconds / bridge_cumulative_seconds : 0.0
    );
    reb_simulation_free(ias15);
    free_bridge_state(&bridge_state);
    return EXIT_SUCCESS;
}
