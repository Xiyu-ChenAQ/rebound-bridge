/* SPDX-License-Identifier: GPL-3.0-or-later */

/**
 * compare_solar_system_ias15.c - Bridge vs IAS15 validation for the solar-system demo.
 *
 * Produces a CSV with bridge and IAS15 diagnostics sampled at identical times.
 */

#include <math.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

#include "rebound_bridge.h"

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

static struct reb_particle particle_from_state(
    double m,
    double x,
    double y,
    double z,
    double vx,
    double vy,
    double vz
) {
    struct reb_particle p;
    memset(&p, 0, sizeof(p));
    p.m = m;
    p.x = x;
    p.y = y;
    p.z = z;
    p.vx = vx;
    p.vy = vy;
    p.vz = vz;
    return p;
}

static struct reb_particle particle_add(const struct reb_particle* host, const struct reb_particle* rel) {
    struct reb_particle p = *rel;
    p.x += host->x;
    p.y += host->y;
    p.z += host->z;
    p.vx += host->vx;
    p.vy += host->vy;
    p.vz += host->vz;
    return p;
}

static int ensure_directory(const char* path) {
#ifdef _WIN32
    if (_mkdir(path) == 0) return 0;
#else
    if (mkdir(path, 0777) == 0) return 0;
#endif
    if (errno == EEXIST) return 0;
    return -1;
}

static int ensure_output_tree(void) {
    if (ensure_directory("validation") != 0) return -1;
    if (ensure_directory("validation/solar_system_ias15") != 0) return -1;
    if (ensure_directory("validation/solar_system_ias15/out") != 0) return -1;
    return 0;
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

static void populate_physical_initial_bodies(struct reb_particle bodies[BODY_COUNT]) {
    bodies[BODY_SUN] = particle_from_state(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    bodies[BODY_MERCURY] = particle_from_state(
        1.660120e-7,
        -2.15201304377667713e-01, -4.09207615700272087e-01, -1.37032695919236600e-02,
        7.02443404629697010e+00, -4.28716799788450942e+00, -9.94632701667513053e-01
    );
    bodies[BODY_VENUS] = particle_from_state(
        2.4478383e-6,
        8.88772496163565545e-02, -7.21762382317080675e-01, -1.50440565750619305e-02,
        7.28255910865529810e+00, 8.76733642322772511e-01, -4.08156573658943866e-01
    );
    bodies[BODY_EARTH] = particle_from_state(
        3.0034896e-6,
        -1.74269758548305093e-01, 9.67785674324191270e-01, -5.68660813624550973e-05,
        -6.28653514426191418e+00, -1.13726197010495311e+00, 1.00945156455726687e-04
    );
    bodies[BODY_MOON] = particle_from_state(
        3.6940e-8,
        -1.74269758548305093e-01 + 9.53035527126304600e-04,
        9.67785674324191270e-01 + 2.17512105974572998e-03,
        -5.68660813624550973e-05 + 2.09676507944608496e-04,
        -6.28653514426191418e+00 - 2.09285374018137460e-01,
        -1.13726197010495311e+00 + 8.77025716574588027e-02,
        1.00945156455726687e-04 + 1.16388188872971468e-03
    );
    bodies[BODY_MARS] = particle_from_state(
        3.227151e-7,
        3.40579676862215075e-01, -1.38700201594525407e+00, -3.74172267877010828e-02,
        5.15711521620981550e+00, 1.65830921706158052e+00, -9.17089273335817273e-02
    );
    bodies[BODY_JUPITER] = particle_from_state(
        9.547919e-4,
        -1.69400303018584597e+00, 4.92888271463623262e+00, 1.74262399033241998e-02,
        -2.64098245722763814e+00, -7.67587434003239921e-01, 6.22754635347801122e-02
    );
    bodies[BODY_IO] = particle_from_state(
        2.388e-8,
        -1.69400303018584597e+00 + 2.48120485765584205e-03,
        4.92888271463623262e+00 - 1.34556438159649208e-03,
        1.74262399033241998e-02 - 1.31630007347439604e-05,
        -2.64098245722763814e+00 + 1.75341752994026812e+00,
        -7.67587434003239921e-01 + 3.20163492363371871e+00,
        6.22754635347801122e-02 + 1.39542430489110164e-01
    );
    bodies[BODY_EUROPA] = particle_from_state(
        2.527e-8,
        -1.69400303018584597e+00 + 5.51813635549570501e-04,
        4.92888271463623262e+00 - 4.47551321186325738e-03,
        1.74262399033241998e-02 - 1.27610904592036212e-04,
        -2.64098245722763814e+00 + 2.85755306406471643e+00,
        -7.67587434003239921e-01 + 3.72439324002732552e-01,
        6.22754635347801122e-02 + 7.16871253077594472e-02
    );
    bodies[BODY_GANYMEDE] = particle_from_state(
        7.804e-8,
        -1.69400303018584597e+00 + 6.76397567656043287e-03,
        4.92888271463623262e+00 - 2.29852901298470890e-03,
        1.74262399033241998e-02 + 1.01873701583559195e-05,
        -2.64098245722763814e+00 + 7.40898140660848692e-01,
        -7.67587434003239921e-01 + 2.17447712709783314e+00,
        6.22754635347801122e-02 + 9.37799108604124437e-02
    );
    bodies[BODY_SATURN] = particle_from_state(
        2.8588567e-4,
        9.50734336522187107e+00, 2.57738457802152210e-01, -3.82935528507753087e-01,
        -1.69141631745304527e-01, 2.03249650346815258e+00, -2.86682018818476998e-02
    );
    bodies[BODY_URANUS] = particle_from_state(
        4.366244e-5,
        9.88031797519427535e+00, 1.68000154098179486e+01, -6.57197991373143220e-02,
        -1.25149298434788658e+00, 6.61240423623731699e-01, 1.86840216329132620e-02
    );
    bodies[BODY_NEPTUNE] = particle_from_state(
        5.151389e-5,
        2.98721200408889800e+01, 5.18939463241125942e-01, -6.99067175203629043e-01,
        -2.98764379233138963e-02, 1.15313346010349438e+00, -2.30704057616862752e-02
    );
}

static struct reb_simulation* make_main_solar_system(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    reb_simulation_add_fmt(sim, "m", 1.0);
    add_cartesian(sim, 1.660120e-7, -2.15201304377667713e-01, -4.09207615700272087e-01, -1.37032695919236600e-02, 7.02443404629697010e+00, -4.28716799788450942e+00, -9.94632701667513053e-01);
    add_cartesian(sim, 2.4478383e-6, 8.88772496163565545e-02, -7.21762382317080675e-01, -1.50440565750619305e-02, 7.28255910865529810e+00, 8.76733642322772511e-01, -4.08156573658943866e-01);
    add_cartesian(sim, 3.0034896e-6 + 3.6940e-8, -1.74269758548305093e-01, 9.67785674324191270e-01, -5.68660813624550973e-05, -6.28653514426191418e+00, -1.13726197010495311e+00, 1.00945156455726687e-04);
    add_cartesian(sim, 3.227151e-7, 3.40579676862215075e-01, -1.38700201594525407e+00, -3.74172267877010828e-02, 5.15711521620981550e+00, 1.65830921706158052e+00, -9.17089273335817273e-02);
    add_cartesian(sim, 9.547919e-4 + 2.388e-8 + 2.527e-8 + 7.804e-8, -1.69400303018584597e+00, 4.92888271463623262e+00, 1.74262399033241998e-02, -2.64098245722763814e+00, -7.67587434003239921e-01, 6.22754635347801122e-02);
    add_cartesian(sim, 2.8588567e-4, 9.50734336522187107e+00, 2.57738457802152210e-01, -3.82935528507753087e-01, -1.69141631745304527e-01, 2.03249650346815258e+00, -2.86682018818476998e-02);
    add_cartesian(sim, 4.366244e-5, 9.88031797519427535e+00, 1.68000154098179486e+01, -6.57197991373143220e-02, -1.25149298434788658e+00, 6.61240423623731699e-01, 1.86840216329132620e-02);
    add_cartesian(sim, 5.151389e-5, 2.98721200408889800e+01, 5.18939463241125942e-01, -6.99067175203629043e-01, -2.98764379233138963e-02, 1.15313346010349438e+00, -2.30704057616862752e-02);

    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_sub_earth_moon_horizons(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    add_cartesian(sim, 3.0034896e-6, -1.17223722782580293e-05, -2.67540695880473792e-05, -2.57902881285356581e-06, 2.57421785107276755e-03, -1.07874487935336202e-03, -1.43157903344392045e-05);
    add_cartesian(sim, 3.6940e-8, 9.53035527126304600e-04, 2.17512105974572998e-03, 2.09676507944608496e-04, -2.09285374018137460e-01, 8.77025716574588027e-02, 1.16388188872971468e-03);

    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_sub_jovian_horizons(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;

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

static struct bridge_state make_bridge_state(double dt_outer, double dt_earth_moon, double dt_jovian) {
    struct bridge_state state;
    memset(&state, 0, sizeof(state));

    state.main_sim = make_main_solar_system();
    state.earth_moon = make_sub_earth_moon_horizons();
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
    if (!sim) return NULL;

    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_ias15(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;

    struct reb_particle bodies[BODY_COUNT];
    populate_physical_initial_bodies(bodies);
    for (int i = 0; i < BODY_COUNT; i++) {
        reb_simulation_add(sim, bodies[i]);
    }

    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_particle barycenter_of_indices(const struct reb_particle* bodies, const int* indices, int count) {
    struct reb_particle com;
    memset(&com, 0, sizeof(com));

    double mt = 0.0;
    for (int i = 0; i < count; i++) {
        const struct reb_particle* p = &bodies[indices[i]];
        mt += p->m;
        com.x += p->m * p->x;
        com.y += p->m * p->y;
        com.z += p->m * p->z;
        com.vx += p->m * p->vx;
        com.vy += p->m * p->vy;
        com.vz += p->m * p->vz;
    }

    if (mt > 0.0) {
        com.x /= mt;
        com.y /= mt;
        com.z /= mt;
        com.vx /= mt;
        com.vy /= mt;
        com.vz /= mt;
    }
    com.m = mt;
    return com;
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
    bodies[BODY_EARTH] = particle_add(host_em, &earth_moon->particles[EM_EARTH]);
    bodies[BODY_MOON] = particle_add(host_em, &earth_moon->particles[EM_MOON]);
    bodies[BODY_MARS] = main_sim->particles[MAIN_MARS];
    bodies[BODY_JUPITER] = particle_add(host_jov, &jovian->particles[JOVIAN_JUPITER]);
    bodies[BODY_IO] = particle_add(host_jov, &jovian->particles[JOVIAN_IO]);
    bodies[BODY_EUROPA] = particle_add(host_jov, &jovian->particles[JOVIAN_EUROPA]);
    bodies[BODY_GANYMEDE] = particle_add(host_jov, &jovian->particles[JOVIAN_GANYMEDE]);
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
        metrics.total_energy += 0.5 * p->m * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
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
    metrics.earth_moon_barycenter = barycenter_of_indices(bodies, em_ids, 2);
    metrics.jovian_barycenter = barycenter_of_indices(bodies, jov_ids, 4);
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

static void write_csv_header(FILE* fp) {
    fprintf(fp,
        "time_yr,"
        "bridge_energy_rel,ias15_energy_rel,energy_rel_diff,"
        "bridge_l_rel,ias15_l_rel,l_rel_diff,"
        "bridge_em_distance_au,ias15_em_distance_au,em_distance_diff_au,"
        "bridge_laplace_deg,ias15_laplace_deg,laplace_diff_deg,"
        "earth_barycenter_error_au,jupiter_barycenter_error_au,"
        "earth_phase_diff_deg,jupiter_phase_diff_deg,"
        "earth_radius_diff_au,jupiter_radius_diff_au\n"
    );
}

static void write_csv_row(
    FILE* fp,
    double time_yr,
    double bridge_energy_rel,
    double ias15_energy_rel,
    double bridge_l_rel,
    double ias15_l_rel,
    double bridge_em_distance,
    double ias15_em_distance,
    double bridge_laplace_deg,
    double ias15_laplace_deg,
    double earth_barycenter_error_au,
    double jupiter_barycenter_error_au,
    double earth_phase_diff_deg,
    double jupiter_phase_diff_deg,
    double earth_radius_diff_au,
    double jupiter_radius_diff_au
) {
    fprintf(fp,
        "%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
        time_yr,
        bridge_energy_rel,
        ias15_energy_rel,
        bridge_energy_rel - ias15_energy_rel,
        bridge_l_rel,
        ias15_l_rel,
        bridge_l_rel - ias15_l_rel,
        bridge_em_distance,
        ias15_em_distance,
        bridge_em_distance - ias15_em_distance,
        bridge_laplace_deg,
        ias15_laplace_deg,
        bridge_laplace_deg - ias15_laplace_deg,
        earth_barycenter_error_au,
        jupiter_barycenter_error_au,
        earth_phase_diff_deg,
        jupiter_phase_diff_deg,
        earth_radius_diff_au,
        jupiter_radius_diff_au
    );
}

static int parse_double_value(const char* text, double* out) {
    char* end = NULL;
    errno = 0;
    const double value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0') return -1;
    *out = value;
    return 0;
}

static int parse_int_value(const char* text, int* out) {
    char* end = NULL;
    errno = 0;
    const long value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value <= 0 || value > INT_MAX) return -1;
    *out = (int)value;
    return 0;
}

static void print_usage(const char* argv0) {
    fprintf(stderr, "usage: %s [--years N] [--samples N] [--dt-outer-days N]\n", argv0);
}

int main(int argc, char** argv) {
    double years = 2000.0;
    int samples = 2000;
    double dt_outer_days = 1.0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--years") == 0) {
            if (i + 1 >= argc || parse_double_value(argv[++i], &years) != 0 || years <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--samples") == 0) {
            if (i + 1 >= argc || parse_int_value(argv[++i], &samples) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--dt-outer-days") == 0) {
            if (i + 1 >= argc || parse_double_value(argv[++i], &dt_outer_days) != 0 || dt_outer_days <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    const double dt_outer = dt_outer_days / 365.25;
    const double dt_earth_moon = dt_outer / 20.0;
    const double dt_jovian = dt_outer / 50.0;
    const char* out_path = "validation/solar_system_ias15/out/solar_system_ias15_compare.csv";

    struct bridge_state bridge_state = make_bridge_state(dt_outer, dt_earth_moon, dt_jovian);
    if (!bridge_state.bridge || !bridge_state.main_sim || !bridge_state.earth_moon || !bridge_state.jovian) {
        fprintf(stderr, "failed to create bridge state\n");
        free_bridge_state(&bridge_state);
        return EXIT_FAILURE;
    }

    struct reb_simulation* ias15 = make_ias15_reference();
    if (!ias15) {
        fprintf(stderr, "failed to create IAS15 reference state\n");
        free_bridge_state(&bridge_state);
        return EXIT_FAILURE;
    }

    if (ensure_output_tree() != 0) {
        fprintf(stderr, "failed to create output directories\n");
        reb_simulation_free(ias15);
        free_bridge_state(&bridge_state);
        return EXIT_FAILURE;
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
    const struct reb_vec3d bridge_l0 = bridge0.total_l;
    const struct reb_vec3d ias15_l0 = ias150.total_l;

    for (int i = 0; i < samples; i++) {
        const double target = years * (double)(i + 1) / (double)samples;

        if (reb_bridge_advance(bridge_state.bridge, target - bridge_state.main_sim->t) != 0) {
            fprintf(stderr, "bridge advance failed at sample %d\n", i + 1);
            fclose(fp);
            reb_simulation_free(ias15);
            free_bridge_state(&bridge_state);
            return EXIT_FAILURE;
        }

        if (reb_simulation_integrate(ias15, target) != REB_STATUS_SUCCESS) {
            fprintf(stderr, "IAS15 integrate failed at sample %d\n", i + 1);
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
        const double bridge_em_distance = distance_particles(&bridge_bodies[BODY_EARTH], &bridge_bodies[BODY_MOON]);
        const double ias15_em_distance = distance_particles(&ias15_bodies[BODY_EARTH], &ias15_bodies[BODY_MOON]);
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
            target,
            bridge_energy_rel,
            ias15_energy_rel,
            bridge_l_rel,
            ias15_l_rel,
            bridge_em_distance,
            ias15_em_distance,
            bridge_laplace_deg,
            ias15_laplace_deg,
            earth_barycenter_error_au,
            jupiter_barycenter_error_au,
            earth_phase_diff_deg,
            jupiter_phase_diff_deg,
            earth_radius_diff_au,
            jupiter_radius_diff_au
        );
    }

    fclose(fp);
    reb_simulation_free(ias15);
    free_bridge_state(&bridge_state);
    return EXIT_SUCCESS;
}
