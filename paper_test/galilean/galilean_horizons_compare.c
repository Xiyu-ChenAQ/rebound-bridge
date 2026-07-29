/* Local real Galilean-satellite Bridge-vs-IAS15 pilot using JPL Horizons states. */

#include <errno.h>
#include <math.h>
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
    BODY_SUN = 0,
    BODY_JUPITER = 1,
    BODY_IO = 2,
    BODY_EUROPA = 3,
    BODY_GANYMEDE = 4,
    BODY_CALLISTO = 5,
    BODY_COUNT = 6
};

struct run_options {
    double years;
    int samples;
    double dt_outer_days;
    int inner_ratio;
    enum reb_bridge_method method;
    const char* method_name;
    const char* csv_path;
};

struct bridge_state {
    struct reb_simulation* main_sim;
    struct reb_simulation* sub_sim;
    struct reb_bridge* bridge;
};

struct system_metrics {
    int ok;
    int escaped;
    double moon_error[4];
    double max_moon_error;
    double laplace_deg;
    double period_ratio_io_europa;
    double period_ratio_europa_ganymede;
    double max_eccentricity;
    double max_distance_au;
};

static const double R_JUPITER_AU = 0.000477894503;

static double norm3(double x, double y, double z) {
    return sqrt(x * x + y * y + z * z);
}

static double distance_particles(const struct reb_particle* a, const struct reb_particle* b) {
    return norm3(a->x - b->x, a->y - b->y, a->z - b->z);
}

static double angle_diff_rad(double a, double b) {
    double d = a - b;
    while (d > M_PI) d -= 2.0 * M_PI;
    while (d < -M_PI) d += 2.0 * M_PI;
    return d;
}

static double wrap_rad(double a) {
    while (a > M_PI) a -= 2.0 * M_PI;
    while (a < -M_PI) a += 2.0 * M_PI;
    return a;
}

static double relative_position_delta(
    const struct reb_particle* a_primary,
    const struct reb_particle* a_secondary,
    const struct reb_particle* b_primary,
    const struct reb_particle* b_secondary
) {
    return norm3(
        (a_secondary->x - a_primary->x) - (b_secondary->x - b_primary->x),
        (a_secondary->y - a_primary->y) - (b_secondary->y - b_primary->y),
        (a_secondary->z - a_primary->z) - (b_secondary->z - b_primary->z)
    );
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
        "usage: %s [--method order2|eos_lf4_2] [--years N] [--samples N] [--dt-outer-days N] [--inner-ratio N] --csv PATH\n",
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

static struct reb_particle galilean_host_com(void) {
    struct reb_particle com;
    memset(&com, 0, sizeof(com));
    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        const struct reb_particle p = particle_from_state(&GALILEAN_BODIES[i]);
        com.m += p.m;
        com.x += p.m * p.x;
        com.y += p.m * p.y;
        com.z += p.m * p.z;
        com.vx += p.m * p.vx;
        com.vy += p.m * p.vy;
        com.vz += p.m * p.vz;
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

static struct reb_simulation* make_main_sim(double dt_outer) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    sim->dt = dt_outer;

    struct reb_particle sun;
    memset(&sun, 0, sizeof(sun));
    sun.m = 1.0;
    const struct reb_particle host = galilean_host_com();
    add_cartesian(sim, &sun);
    add_cartesian(sim, &host);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_sub_sim(double dt_inner) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    sim->dt = dt_inner;

    const struct reb_particle host = galilean_host_com();
    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        const struct reb_particle absolute = particle_from_state(&GALILEAN_BODIES[i]);
        const struct reb_particle relative = subtract_particle(&absolute, &host);
        add_cartesian(sim, &relative);
    }
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct bridge_state make_bridge_state(const struct run_options* options) {
    struct bridge_state state;
    memset(&state, 0, sizeof(state));
    const double dt_outer = options->dt_outer_days / 365.25;
    const double dt_inner = dt_outer / (double)options->inner_ratio;

    state.main_sim = make_main_sim(dt_outer);
    state.sub_sim = make_sub_sim(dt_inner);
    if (!state.main_sim || !state.sub_sim) return state;

    state.bridge = reb_bridge_create(state.main_sim, dt_outer);
    if (!state.bridge) return state;
    reb_bridge_set_owns_main(state.bridge, 1);
    if (reb_bridge_set_method(state.bridge, options->method) != 0) return state;
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

static struct reb_simulation* make_full_ias15(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_ias15(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;

    struct reb_particle sun;
    memset(&sun, 0, sizeof(sun));
    sun.m = 1.0;
    add_cartesian(sim, &sun);
    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        const struct reb_particle p = particle_from_state(&GALILEAN_BODIES[i]);
        add_cartesian(sim, &p);
    }
    reb_simulation_move_to_com(sim);
    return sim;
}

static void fill_bridge_bodies(const struct bridge_state* state, struct reb_particle bodies[BODY_COUNT]) {
    const struct reb_particle host = state->main_sim->particles[1];
    bodies[BODY_SUN] = state->main_sim->particles[0];
    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        bodies[i + 1] = add_particle(&state->sub_sim->particles[i], &host);
    }
}

static void fill_full_bodies(const struct reb_simulation* sim, struct reb_particle bodies[BODY_COUNT]) {
    for (int i = 0; i < BODY_COUNT; i++) bodies[i] = sim->particles[i];
}

static double total_energy(const struct reb_particle bodies[BODY_COUNT], double G) {
    double energy = 0.0;
    for (int i = 0; i < BODY_COUNT; i++) {
        const double v2 = bodies[i].vx * bodies[i].vx
            + bodies[i].vy * bodies[i].vy
            + bodies[i].vz * bodies[i].vz;
        energy += 0.5 * bodies[i].m * v2;
    }
    for (int i = 0; i < BODY_COUNT; i++) {
        for (int j = i + 1; j < BODY_COUNT; j++) {
            const double r = distance_particles(&bodies[i], &bodies[j]);
            energy -= G * bodies[i].m * bodies[j].m / r;
        }
    }
    return energy;
}

static struct system_metrics compute_metrics(
    const struct reb_particle bodies[BODY_COUNT],
    const struct reb_particle* reference_bodies,
    double G
) {
    struct system_metrics out;
    memset(&out, 0, sizeof(out));
    out.ok = 1;

    struct reb_orbit orbits[4];
    for (int i = 0; i < 4; i++) {
        int err = 0;
        const int body_id = BODY_IO + i;
        orbits[i] = reb_orbit_from_particle_err(G, bodies[body_id], bodies[BODY_JUPITER], &err);
        if (err || !isfinite(orbits[i].a) || !isfinite(orbits[i].e) || !isfinite(orbits[i].l)) {
            out.ok = 0;
            return out;
        }
        const double distance = distance_particles(&bodies[BODY_JUPITER], &bodies[body_id]);
        if (distance > out.max_distance_au) out.max_distance_au = distance;
        if (orbits[i].e > out.max_eccentricity) out.max_eccentricity = orbits[i].e;
        if (distance < R_JUPITER_AU || orbits[i].e >= 1.0) out.escaped = 1;
        if (reference_bodies) {
            out.moon_error[i] = relative_position_delta(
                &bodies[BODY_JUPITER],
                &bodies[body_id],
                &reference_bodies[BODY_JUPITER],
                &reference_bodies[body_id]
            );
            if (out.moon_error[i] > out.max_moon_error) out.max_moon_error = out.moon_error[i];
        }
    }

    out.period_ratio_io_europa = orbits[1].P / orbits[0].P;
    out.period_ratio_europa_ganymede = orbits[2].P / orbits[1].P;
    out.laplace_deg = wrap_rad(orbits[0].l - 3.0 * orbits[1].l + 2.0 * orbits[2].l) * 180.0 / M_PI;
    return out;
}

static int parse_args(int argc, char** argv, struct run_options* options) {
    options->years = 20.0;
    options->samples = 200;
    options->dt_outer_days = 0.5;
    options->inner_ratio = 40;
    options->method = REB_BRIDGE_METHOD_EOS_LF4_2;
    options->method_name = "eos_lf4_2";
    options->csv_path = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--method") == 0) {
            if (++i >= argc) return -1;
            if (strcmp(argv[i], "order2") == 0) {
                options->method = REB_BRIDGE_METHOD_ORDER2;
                options->method_name = "order2";
            } else if (strcmp(argv[i], "eos_lf4_2") == 0) {
                options->method = REB_BRIDGE_METHOD_EOS_LF4_2;
                options->method_name = "eos_lf4_2";
            } else {
                return -1;
            }
        } else if (strcmp(argv[i], "--years") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->years) != 0 || options->years <= 0.0) return -1;
        } else if (strcmp(argv[i], "--samples") == 0) {
            if (++i >= argc || parse_int_arg(argv[i], &options->samples) != 0) return -1;
        } else if (strcmp(argv[i], "--dt-outer-days") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->dt_outer_days) != 0 || options->dt_outer_days <= 0.0) return -1;
        } else if (strcmp(argv[i], "--inner-ratio") == 0) {
            if (++i >= argc || parse_int_arg(argv[i], &options->inner_ratio) != 0) return -1;
        } else if (strcmp(argv[i], "--csv") == 0) {
            if (++i >= argc) return -1;
            options->csv_path = argv[i];
        } else {
            return -1;
        }
    }

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

    struct bridge_state bridge = make_bridge_state(&options);
    struct reb_simulation* ias15 = make_full_ias15();
    if (!bridge.bridge || !bridge.main_sim || !bridge.sub_sim || !ias15) {
        free_bridge_state(&bridge);
        if (ias15) reb_simulation_free(ias15);
        fclose(csv);
        fprintf(stderr, "failed to initialize simulations\n");
        return 1;
    }

    fprintf(
        csv,
        "time_yr,bridge_step_seconds,ias15_step_seconds,bridge_cumulative_seconds,ias15_cumulative_seconds,cumulative_speedup,"
        "bridge_ok,ias15_ok,bridge_escaped,ias15_escaped,"
        "io_error_au,europa_error_au,ganymede_error_au,callisto_error_au,max_moon_error_au,"
        "bridge_laplace_deg,ias15_laplace_deg,laplace_diff_deg,"
        "bridge_period_ratio_io_europa,ias15_period_ratio_io_europa,period_ratio_io_europa_diff,"
        "bridge_period_ratio_europa_ganymede,ias15_period_ratio_europa_ganymede,period_ratio_europa_ganymede_diff,"
        "bridge_max_e,ias15_max_e,bridge_max_distance_au,ias15_max_distance_au,"
        "bridge_energy,ias15_energy,bridge_energy_rel,ias15_energy_rel,energy_rel_diff\n"
    );

    printf("Galilean Horizons Bridge-vs-IAS15 pilot\n");
    printf("epoch=%s source=%s\n", GALILEAN_HORIZONS_EPOCH, GALILEAN_HORIZONS_SOURCE);
    printf("method=%s years=%.17g samples=%d dt_outer_days=%.17g inner_dt_days=%.17g\n",
        options.method_name,
        options.years,
        options.samples,
        options.dt_outer_days,
        options.dt_outer_days / (double)options.inner_ratio
    );

    double bridge_total_seconds = 0.0;
    double ias15_total_seconds = 0.0;
    double max_moon_error = 0.0;
    double max_laplace_diff_deg = 0.0;
    double max_period_ratio_io_europa_diff = 0.0;
    double max_period_ratio_europa_ganymede_diff = 0.0;
    int bridge_escaped = 0;
    int ias15_escaped = 0;

    struct reb_particle initial_bridge_bodies[BODY_COUNT];
    struct reb_particle initial_ias15_bodies[BODY_COUNT];
    fill_bridge_bodies(&bridge, initial_bridge_bodies);
    fill_full_bodies(ias15, initial_ias15_bodies);
    const double initial_bridge_energy = total_energy(initial_bridge_bodies, bridge.main_sim->G);
    const double initial_ias15_energy = total_energy(initial_ias15_bodies, ias15->G);
    double max_abs_bridge_energy_rel = 0.0;
    double max_abs_ias15_energy_rel = 0.0;

    for (int sample = 1; sample <= options.samples; sample++) {
        const double target = options.years * (double)sample / (double)options.samples;

        const double bridge_start = monotonic_seconds();
        if (reb_bridge_advance(bridge.bridge, target - bridge.main_sim->t) != 0) {
            fprintf(stderr, "bridge integration failed at %.17g yr\n", target);
            free_bridge_state(&bridge);
            reb_simulation_free(ias15);
            fclose(csv);
            return 1;
        }
        const double bridge_step_seconds = monotonic_seconds() - bridge_start;
        bridge_total_seconds += bridge_step_seconds;

        const double ias15_start = monotonic_seconds();
        if (reb_simulation_integrate(ias15, target) != REB_STATUS_SUCCESS) {
            fprintf(stderr, "ias15 integration failed at %.17g yr\n", target);
            free_bridge_state(&bridge);
            reb_simulation_free(ias15);
            fclose(csv);
            return 1;
        }
        const double ias15_step_seconds = monotonic_seconds() - ias15_start;
        ias15_total_seconds += ias15_step_seconds;

        struct reb_particle bridge_bodies[BODY_COUNT];
        struct reb_particle ias15_bodies[BODY_COUNT];
        fill_bridge_bodies(&bridge, bridge_bodies);
        fill_full_bodies(ias15, ias15_bodies);

        const struct system_metrics ias15_metrics = compute_metrics(ias15_bodies, NULL, ias15->G);
        const struct system_metrics bridge_metrics = compute_metrics(bridge_bodies, ias15_bodies, bridge.main_sim->G);
        const double bridge_energy = total_energy(bridge_bodies, bridge.main_sim->G);
        const double ias15_energy = total_energy(ias15_bodies, ias15->G);
        const double bridge_energy_rel = (bridge_energy - initial_bridge_energy) / fabs(initial_bridge_energy);
        const double ias15_energy_rel = (ias15_energy - initial_ias15_energy) / fabs(initial_ias15_energy);
        const double energy_rel_diff = bridge_energy_rel - ias15_energy_rel;
        const double laplace_diff_deg = fabs(angle_diff_rad(
            bridge_metrics.laplace_deg * M_PI / 180.0,
            ias15_metrics.laplace_deg * M_PI / 180.0
        )) * 180.0 / M_PI;
        const double period_ratio_io_europa_diff =
            fabs(bridge_metrics.period_ratio_io_europa - ias15_metrics.period_ratio_io_europa);
        const double period_ratio_europa_ganymede_diff =
            fabs(bridge_metrics.period_ratio_europa_ganymede - ias15_metrics.period_ratio_europa_ganymede);

        if (bridge_metrics.max_moon_error > max_moon_error) max_moon_error = bridge_metrics.max_moon_error;
        if (laplace_diff_deg > max_laplace_diff_deg) max_laplace_diff_deg = laplace_diff_deg;
        if (period_ratio_io_europa_diff > max_period_ratio_io_europa_diff) {
            max_period_ratio_io_europa_diff = period_ratio_io_europa_diff;
        }
        if (period_ratio_europa_ganymede_diff > max_period_ratio_europa_ganymede_diff) {
            max_period_ratio_europa_ganymede_diff = period_ratio_europa_ganymede_diff;
        }
        if (bridge_metrics.escaped) bridge_escaped = 1;
        if (ias15_metrics.escaped) ias15_escaped = 1;
        if (fabs(bridge_energy_rel) > max_abs_bridge_energy_rel) {
            max_abs_bridge_energy_rel = fabs(bridge_energy_rel);
        }
        if (fabs(ias15_energy_rel) > max_abs_ias15_energy_rel) {
            max_abs_ias15_energy_rel = fabs(ias15_energy_rel);
        }

        fprintf(
            csv,
            "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%d,%d,%d,%d,"
            "%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,%.17g,%.17g\n",
            target,
            bridge_step_seconds,
            ias15_step_seconds,
            bridge_total_seconds,
            ias15_total_seconds,
            bridge_total_seconds > 0.0 ? ias15_total_seconds / bridge_total_seconds : 0.0,
            bridge_metrics.ok,
            ias15_metrics.ok,
            bridge_metrics.escaped,
            ias15_metrics.escaped,
            bridge_metrics.moon_error[0],
            bridge_metrics.moon_error[1],
            bridge_metrics.moon_error[2],
            bridge_metrics.moon_error[3],
            bridge_metrics.max_moon_error,
            bridge_metrics.laplace_deg,
            ias15_metrics.laplace_deg,
            laplace_diff_deg,
            bridge_metrics.period_ratio_io_europa,
            ias15_metrics.period_ratio_io_europa,
            period_ratio_io_europa_diff,
            bridge_metrics.period_ratio_europa_ganymede,
            ias15_metrics.period_ratio_europa_ganymede,
            period_ratio_europa_ganymede_diff,
            bridge_metrics.max_eccentricity,
            ias15_metrics.max_eccentricity,
            bridge_metrics.max_distance_au,
            ias15_metrics.max_distance_au,
            bridge_energy,
            ias15_energy,
            bridge_energy_rel,
            ias15_energy_rel,
            energy_rel_diff
        );
        fflush(csv);
    }

    printf("bridge_total_seconds=%.9f\n", bridge_total_seconds);
    printf("ias15_total_seconds=%.9f\n", ias15_total_seconds);
    printf("ias15_over_bridge_speedup=%.9f\n", bridge_total_seconds > 0.0 ? ias15_total_seconds / bridge_total_seconds : 0.0);
    printf("max_moon_error_au=%.17g\n", max_moon_error);
    printf("max_laplace_diff_deg=%.17g\n", max_laplace_diff_deg);
    printf("max_period_ratio_io_europa_diff=%.17g\n", max_period_ratio_io_europa_diff);
    printf("max_period_ratio_europa_ganymede_diff=%.17g\n", max_period_ratio_europa_ganymede_diff);
    printf("max_abs_bridge_energy_rel=%.17g\n", max_abs_bridge_energy_rel);
    printf("max_abs_ias15_energy_rel=%.17g\n", max_abs_ias15_energy_rel);
    printf("bridge_escaped=%d\n", bridge_escaped);
    printf("ias15_escaped=%d\n", ias15_escaped);
    printf("wrote %s\n", options.csv_path);

    free_bridge_state(&bridge);
    reb_simulation_free(ias15);
    fclose(csv);
    return 0;
}
