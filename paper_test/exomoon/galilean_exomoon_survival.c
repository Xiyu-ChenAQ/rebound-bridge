/* Galilean-like resonant moon chain survival around a Sun-like host star. */

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "galilean_horizons_state.h"
#include "rebound_bridge.h"
#include "validation_timer.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum full_body {
    BODY_STAR = 0,
    BODY_JUPITER = 1,
    BODY_IO = 2,
    BODY_EUROPA = 3,
    BODY_GANYMEDE = 4,
    BODY_CALLISTO = 5,
    BODY_COUNT = 6
};

struct run_options {
    const char* integrator_name;
    double a_planet_au;
    double planet_eccentricity;
    double planet_phase_deg;
    double years;
    double dt_outer_days;
    int inner_ratio;
    int samples;
    const char* csv_path;
};

struct bridge_state {
    struct reb_simulation* main_sim;
    struct reb_simulation* sub_sim;
    struct reb_bridge* bridge;
};

struct chain_metrics {
    int ok;
    int escaped;
    int collided;
    double max_eccentricity;
    double max_distance_au;
    double laplace_deg;
    double inner_resonant_angle_deg;
    double outer_resonant_angle_deg;
    double laplace_detuning_rad_per_yr;
    double period_ratio_io_europa;
    double period_ratio_europa_ganymede;
};

static const double R_JUPITER_AU = 0.000477894503;

static double norm3(double x, double y, double z) {
    return sqrt(x * x + y * y + z * z);
}

static double distance_particles(const struct reb_particle* a, const struct reb_particle* b) {
    return norm3(a->x - b->x, a->y - b->y, a->z - b->z);
}

static double wrap_rad(double a) {
    while (a > M_PI) a -= 2.0 * M_PI;
    while (a < -M_PI) a += 2.0 * M_PI;
    return a;
}

static double angle_diff_deg(double a_deg, double b_deg) {
    return fabs(wrap_rad((a_deg - b_deg) * M_PI / 180.0)) * 180.0 / M_PI;
}

static int parse_double_arg(const char* text, double* out) {
    char* end = NULL;
    errno = 0;
    const double value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !isfinite(value)) return -1;
    *out = value;
    return 0;
}

static int parse_int_arg(const char* text, int* out) {
    char* end = NULL;
    errno = 0;
    const long value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value < 1 || value > 1000000L) return -1;
    *out = (int)value;
    return 0;
}

static void print_usage(const char* argv0) {
    fprintf(
        stderr,
        "usage: %s --a-planet-au A --csv PATH [--integrator bridge|ias15] [--planet-e N] [--years N] "
        "[--planet-phase-deg N] [--dt-outer-days N] [--inner-ratio N] "
        "[--samples N]\n",
        argv0
    );
}

static struct reb_particle particle_from_state(const struct galilean_body_state* state) {
    struct reb_particle p;
    memset(&p, 0, sizeof(p));
    p.m = state->mass_msun;
    p.x = state->x_au;
    p.y = state->y_au;
    p.z = state->z_au;
    p.vx = state->vx_au_yr;
    p.vy = state->vy_au_yr;
    p.vz = state->vz_au_yr;
    return p;
}

static void add_cartesian(struct reb_simulation* sim, const struct reb_particle* p) {
    reb_simulation_add_fmt(
        sim,
        "m x y z vx vy vz",
        p->m,
        p->x,
        p->y,
        p->z,
        p->vx,
        p->vy,
        p->vz
    );
}

static struct reb_particle host_com(const struct reb_particle bodies[BODY_COUNT]) {
    struct reb_particle com;
    memset(&com, 0, sizeof(com));
    for (int i = BODY_JUPITER; i < BODY_COUNT; i++) {
        const struct reb_particle* p = &bodies[i];
        com.m += p->m;
        com.x += p->m * p->x;
        com.y += p->m * p->y;
        com.z += p->m * p->z;
        com.vx += p->m * p->vx;
        com.vy += p->m * p->vy;
        com.vz += p->m * p->vz;
    }
    com.x /= com.m;
    com.y /= com.m;
    com.z /= com.m;
    com.vx /= com.m;
    com.vy /= com.m;
    com.vz /= com.m;
    return com;
}

static struct reb_particle subtract_particle(const struct reb_particle* p, const struct reb_particle* origin) {
    struct reb_particle out = *p;
    out.x -= origin->x;
    out.y -= origin->y;
    out.z -= origin->z;
    out.vx -= origin->vx;
    out.vy -= origin->vy;
    out.vz -= origin->vz;
    return out;
}

static struct reb_particle add_particle(const struct reb_particle* p, const struct reb_particle* origin) {
    struct reb_particle out = *p;
    out.x += origin->x;
    out.y += origin->y;
    out.z += origin->z;
    out.vx += origin->vx;
    out.vy += origin->vy;
    out.vz += origin->vz;
    return out;
}

static void rotate_particle_xy(struct reb_particle* p, double angle_rad) {
    const double c = cos(angle_rad);
    const double s = sin(angle_rad);
    const double x = p->x;
    const double y = p->y;
    const double vx = p->vx;
    const double vy = p->vy;
    p->x = c * x - s * y;
    p->y = s * x + c * y;
    p->vx = c * vx - s * vy;
    p->vy = s * vx + c * vy;
}

static void make_horizons_template(struct reb_particle bodies[BODY_COUNT]) {
    memset(bodies, 0, sizeof(struct reb_particle) * BODY_COUNT);
    bodies[BODY_STAR].m = 1.0;
    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        bodies[i + 1] = particle_from_state(&GALILEAN_BODIES[i]);
    }
}

static void make_galilean_like_bodies(
    double a_planet_au,
    double planet_eccentricity,
    double planet_phase_deg,
    struct reb_particle bodies[BODY_COUNT]
) {
    struct reb_particle template_bodies[BODY_COUNT];
    make_horizons_template(template_bodies);
    const struct reb_particle template_host = host_com(template_bodies);

    struct reb_particle target_host;
    memset(&target_host, 0, sizeof(target_host));
    target_host.m = template_host.m;
    target_host.x = a_planet_au * (1.0 - planet_eccentricity);
    target_host.vy = sqrt(
        REB_BRIDGE_G_AU_YR_MSUN *
        (1.0 + target_host.m) *
        (1.0 + planet_eccentricity) /
        (a_planet_au * (1.0 - planet_eccentricity))
    );
    rotate_particle_xy(&target_host, planet_phase_deg * M_PI / 180.0);

    memset(bodies, 0, sizeof(struct reb_particle) * BODY_COUNT);
    bodies[BODY_STAR].m = 1.0;
    for (int i = BODY_JUPITER; i < BODY_COUNT; i++) {
        const struct reb_particle relative = subtract_particle(&template_bodies[i], &template_host);
        bodies[i] = add_particle(&relative, &target_host);
    }
}

static struct chain_metrics compute_chain_metrics(
    const struct reb_particle bodies[BODY_COUNT],
    double G,
    double hill_radius_au
) {
    struct chain_metrics out;
    memset(&out, 0, sizeof(out));
    out.ok = 1;

    struct reb_orbit orbits[4];
    for (int i = 0; i < 4; i++) {
        int err = 0;
        const int body_id = BODY_IO + i;
        const double distance = distance_particles(&bodies[BODY_JUPITER], &bodies[body_id]);
        if (distance > out.max_distance_au) out.max_distance_au = distance;
        if (distance < R_JUPITER_AU) out.collided = 1;
        if (distance > 0.5 * hill_radius_au) out.escaped = 1;

        orbits[i] = reb_orbit_from_particle_err(G, bodies[body_id], bodies[BODY_JUPITER], &err);
        if (err || !isfinite(orbits[i].e) || !isfinite(orbits[i].l)) {
            out.ok = 0;
            out.escaped = 1;
            continue;
        }
        if (orbits[i].e > out.max_eccentricity) out.max_eccentricity = orbits[i].e;
        if (orbits[i].e >= 1.0) out.escaped = 1;
    }

    if (out.ok) {
        out.laplace_deg = wrap_rad(orbits[0].l - 3.0 * orbits[1].l + 2.0 * orbits[2].l) * 180.0 / M_PI;
        /* The Laplace angle is a derived diagnostic; classify the chain using
         * the two constituent 2:1 angles instead. */
        out.inner_resonant_angle_deg = wrap_rad(
            orbits[0].l - 2.0 * orbits[1].l + orbits[1].pomega
        ) * 180.0 / M_PI;
        out.outer_resonant_angle_deg = wrap_rad(
            orbits[1].l - 2.0 * orbits[2].l + orbits[1].pomega
        ) * 180.0 / M_PI;
        const double n_io = 2.0 * M_PI / orbits[0].P;
        const double n_europa = 2.0 * M_PI / orbits[1].P;
        const double n_ganymede = 2.0 * M_PI / orbits[2].P;
        out.laplace_detuning_rad_per_yr = n_io - 3.0 * n_europa + 2.0 * n_ganymede;
        out.period_ratio_io_europa = orbits[1].P / orbits[0].P;
        out.period_ratio_europa_ganymede = orbits[2].P / orbits[1].P;
    }
    return out;
}

static void fill_bridge_bodies(const struct bridge_state* state, struct reb_particle bodies[BODY_COUNT]) {
    const struct reb_particle host = state->main_sim->particles[1];
    bodies[BODY_STAR] = state->main_sim->particles[0];
    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        bodies[i + 1] = add_particle(&state->sub_sim->particles[i], &host);
    }
}

static struct reb_simulation* make_direct_state(const struct reb_particle bodies[BODY_COUNT]) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_ias15(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    for (int i = 0; i < BODY_COUNT; i++) add_cartesian(sim, &bodies[i]);
    reb_simulation_move_to_com(sim);
    return sim;
}

static void fill_direct_bodies(const struct reb_simulation* sim, struct reb_particle bodies[BODY_COUNT]) {
    for (int i = 0; i < BODY_COUNT; i++) bodies[i] = sim->particles[i];
}

static struct bridge_state make_bridge_state(
    const struct reb_particle bodies[BODY_COUNT],
    double dt_outer_days,
    int inner_ratio
) {
    struct bridge_state state;
    memset(&state, 0, sizeof(state));
    const double dt_outer = dt_outer_days / 365.25;
    const double dt_inner = dt_outer / (double)inner_ratio;
    const struct reb_particle host = host_com(bodies);

    state.main_sim = reb_simulation_create();
    state.sub_sim = reb_simulation_create();
    if (!state.main_sim || !state.sub_sim) return state;
    state.main_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    state.sub_sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(state.main_sim) != 0) return state;
    if (reb_bridge_set_integrator_whfast(state.sub_sim) != 0) return state;
    state.main_sim->exact_finish_time = 1;
    state.sub_sim->exact_finish_time = 1;
    state.main_sim->dt = dt_outer;
    state.sub_sim->dt = dt_inner;

    add_cartesian(state.main_sim, &bodies[BODY_STAR]);
    add_cartesian(state.main_sim, &host);
    reb_simulation_move_to_com(state.main_sim);

    for (int i = BODY_JUPITER; i < BODY_COUNT; i++) {
        const struct reb_particle relative = subtract_particle(&bodies[i], &host);
        add_cartesian(state.sub_sim, &relative);
    }
    reb_simulation_move_to_com(state.sub_sim);

    state.bridge = reb_bridge_create(state.main_sim, dt_outer);
    if (!state.bridge) return state;
    reb_bridge_set_owns_main(state.bridge, 1);
    if (reb_bridge_set_method(state.bridge, REB_BRIDGE_METHOD_EOS_LF4_2) != 0) return state;
    if (reb_bridge_add_subsystem(state.bridge, 1, state.sub_sim, dt_inner, 1) != 0) return state;
    return state;
}

static void free_bridge_state(struct bridge_state* state) {
    if (!state) return;
    if (state->bridge) {
        reb_bridge_free(state->bridge);
    } else {
        if (state->main_sim) reb_simulation_free(state->main_sim);
        if (state->sub_sim) reb_simulation_free(state->sub_sim);
    }
    memset(state, 0, sizeof(*state));
}

static int parse_args(int argc, char** argv, struct run_options* options) {
    options->integrator_name = "bridge";
    options->a_planet_au = 0.0;
    options->planet_eccentricity = 0.0;
    options->planet_phase_deg = 0.0;
    options->years = 1000.0;
    options->dt_outer_days = 0.5;
    options->inner_ratio = 40;
    options->samples = 100;
    options->csv_path = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--integrator") == 0) {
            if (++i >= argc) return -1;
            if (strcmp(argv[i], "bridge") != 0 && strcmp(argv[i], "ias15") != 0) return -1;
            options->integrator_name = argv[i];
        } else if (strcmp(argv[i], "--a-planet-au") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->a_planet_au) != 0) return -1;
        } else if (strcmp(argv[i], "--planet-e") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->planet_eccentricity) != 0) return -1;
        } else if (strcmp(argv[i], "--planet-phase-deg") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->planet_phase_deg) != 0) return -1;
        } else if (strcmp(argv[i], "--years") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->years) != 0) return -1;
        } else if (strcmp(argv[i], "--dt-outer-days") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->dt_outer_days) != 0) return -1;
        } else if (strcmp(argv[i], "--inner-ratio") == 0) {
            if (++i >= argc || parse_int_arg(argv[i], &options->inner_ratio) != 0) return -1;
        } else if (strcmp(argv[i], "--samples") == 0) {
            if (++i >= argc || parse_int_arg(argv[i], &options->samples) != 0) return -1;
        } else if (strcmp(argv[i], "--csv") == 0) {
            if (++i >= argc) return -1;
            options->csv_path = argv[i];
        } else {
            return -1;
        }
    }

    if (options->a_planet_au <= 0.0) return -1;
    if (options->planet_eccentricity < 0.0 || options->planet_eccentricity >= 1.0) return -1;
    if (!isfinite(options->planet_phase_deg)) return -1;
    if (options->years <= 0.0 || options->dt_outer_days <= 0.0) return -1;
    return options->csv_path ? 0 : -1;
}

int main(int argc, char** argv) {
    struct run_options options;
    if (parse_args(argc, argv, &options) != 0) {
        print_usage(argv[0]);
        return 2;
    }

    FILE* csv = fopen(options.csv_path, "w");
    if (!csv) {
        fprintf(stderr, "failed to open CSV path: %s\n", options.csv_path);
        return 1;
    }

    struct reb_particle initial_bodies[BODY_COUNT];
    make_galilean_like_bodies(
        options.a_planet_au,
        options.planet_eccentricity,
        options.planet_phase_deg,
        initial_bodies
    );
    const struct reb_particle initial_host = host_com(initial_bodies);
    const double pericentre_au = options.a_planet_au * (1.0 - options.planet_eccentricity);
    const double hill_radius_au = pericentre_au * cbrt(initial_host.m / (3.0 * initial_bodies[BODY_STAR].m));

    struct bridge_state bridge;
    memset(&bridge, 0, sizeof(bridge));
    struct reb_simulation* direct = NULL;
    if (strcmp(options.integrator_name, "bridge") == 0) {
        bridge = make_bridge_state(initial_bodies, options.dt_outer_days, options.inner_ratio);
    } else {
        direct = make_direct_state(initial_bodies);
    }
    if ((strcmp(options.integrator_name, "bridge") == 0 &&
         (!bridge.bridge || !bridge.main_sim || !bridge.sub_sim)) ||
        (strcmp(options.integrator_name, "ias15") == 0 && !direct)) {
        free_bridge_state(&bridge);
        if (direct) reb_simulation_free(direct);
        fclose(csv);
        fprintf(stderr, "failed to initialize %s survival run\n", options.integrator_name);
        return 1;
    }

    struct reb_particle current_bodies[BODY_COUNT];
    if (direct) fill_direct_bodies(direct, current_bodies);
    else fill_bridge_bodies(&bridge, current_bodies);
    const struct chain_metrics initial_metrics =
        compute_chain_metrics(current_bodies, REB_BRIDGE_G_AU_YR_MSUN, hill_radius_au);
    const double initial_laplace_deg = initial_metrics.laplace_deg;

    const uint64_t total_steps = (uint64_t)llround(options.years * 365.25 / options.dt_outer_days);
    const uint64_t output_interval = total_steps / (uint64_t)options.samples > 0
        ? total_steps / (uint64_t)options.samples
        : 1;
    const double dt_years = options.dt_outer_days / 365.25;
    double max_laplace_delta_deg = 0.0;
    double max_eccentricity_seen = 0.0;
    double max_distance_seen_au = 0.0;
    int ever_escaped = 0;
    int ever_collided = 0;

    fprintf(
        csv,
        "a_planet_au,planet_eccentricity,planet_phase_deg,time_yr,step,wall_seconds,hill_radius_au,"
        "max_distance_au,max_distance_over_hill,max_eccentricity,laplace_deg,laplace_delta_deg,"
        "max_laplace_delta_deg,laplace_detuning_rad_per_yr,period_ratio_io_europa,"
        "period_ratio_europa_ganymede,inner_resonant_angle_deg,outer_resonant_angle_deg,"
        "ever_escaped,ever_collided\n"
    );

    printf("Galilean-like exomoon survival run\n");
    printf("epoch_template=%s source=%s\n", GALILEAN_HORIZONS_EPOCH, GALILEAN_HORIZONS_SOURCE);
    printf("integrator=%s a_planet_au=%.17g planet_e=%.17g planet_phase_deg=%.17g years=%.17g dt_outer_days=%.17g inner_ratio=%d hill_radius_au=%.17g\n",
        options.integrator_name,
        options.a_planet_au,
        options.planet_eccentricity,
        options.planet_phase_deg,
        options.years,
        options.dt_outer_days,
        options.inner_ratio,
        hill_radius_au
    );

    const double start = monotonic_seconds();
    for (uint64_t step = 0; step <= total_steps; step++) {
        if (step > 0) {
            const int status = direct
                ? reb_simulation_integrate(direct, (double)step * dt_years)
                : reb_bridge_step(bridge.bridge, dt_years);
            if (status != REB_STATUS_SUCCESS) {
                free_bridge_state(&bridge);
                if (direct) reb_simulation_free(direct);
                fclose(csv);
                fprintf(stderr, "%s step failed at step %llu\n", options.integrator_name, (unsigned long long)step);
                return 1;
            }
        }

        if (step == 0 || step == total_steps || step % output_interval == 0) {
            if (direct) fill_direct_bodies(direct, current_bodies);
            else fill_bridge_bodies(&bridge, current_bodies);
            const struct chain_metrics metrics =
                compute_chain_metrics(current_bodies, REB_BRIDGE_G_AU_YR_MSUN, hill_radius_au);
            const double laplace_delta_deg = angle_diff_deg(metrics.laplace_deg, initial_laplace_deg);
            if (laplace_delta_deg > max_laplace_delta_deg) max_laplace_delta_deg = laplace_delta_deg;
            if (metrics.max_eccentricity > max_eccentricity_seen) max_eccentricity_seen = metrics.max_eccentricity;
            if (metrics.max_distance_au > max_distance_seen_au) max_distance_seen_au = metrics.max_distance_au;
            if (metrics.escaped) ever_escaped = 1;
            if (metrics.collided) ever_collided = 1;

            fprintf(
                csv,
                "%.17g,%.17g,%.17g,%.17g,%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d\n",
                options.a_planet_au,
                options.planet_eccentricity,
                options.planet_phase_deg,
                (double)step * dt_years,
                (unsigned long long)step,
                monotonic_seconds() - start,
                hill_radius_au,
                max_distance_seen_au,
                hill_radius_au > 0.0 ? max_distance_seen_au / hill_radius_au : 0.0,
                max_eccentricity_seen,
                metrics.laplace_deg,
                laplace_delta_deg,
                max_laplace_delta_deg,
                metrics.laplace_detuning_rad_per_yr,
                metrics.period_ratio_io_europa,
                metrics.period_ratio_europa_ganymede,
                metrics.inner_resonant_angle_deg,
                metrics.outer_resonant_angle_deg,
                ever_escaped,
                ever_collided
            );
            fflush(csv);
        }
    }

    printf("wall_seconds=%.9f\n", monotonic_seconds() - start);
    printf("max_distance_over_hill=%.17g\n", hill_radius_au > 0.0 ? max_distance_seen_au / hill_radius_au : 0.0);
    printf("max_eccentricity=%.17g\n", max_eccentricity_seen);
    printf("max_laplace_delta_deg=%.17g\n", max_laplace_delta_deg);
    printf("ever_escaped=%d\n", ever_escaped);
    printf("ever_collided=%d\n", ever_collided);
    printf("wrote %s\n", options.csv_path);

    free_bridge_state(&bridge);
    if (direct) reb_simulation_free(direct);
    fclose(csv);
    return 0;
}
