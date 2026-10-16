/* Local naive global WHFast/SABA vs IAS15 Galilean pilot using JPL Horizons states. */

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
    double dt_days;
    int whfast_corrector;
    int whfast_corrector2;
    const char* integrator_name;
    const char* csv_path;
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

static int parse_nonnegative_int_arg(const char* text, int* out) {
    char* end = NULL;
    errno = 0;
    const long value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value < 0 || value > 1000000L) return -1;
    *out = (int)value;
    return 0;
}

static int valid_whfast_corrector(int corrector) {
    return (
        corrector == 0 ||
        corrector == 3 ||
        corrector == 5 ||
        corrector == 7 ||
        corrector == 11 ||
        corrector == 17
    );
}

static void print_usage(const char* argv0) {
    fprintf(
        stderr,
        "usage: %s [--integrator whfast|saba] [--years N] [--samples N] [--dt-days N] "
        "[--whfast-corrector 0|3|5|7|11|17] [--whfast-corrector2 0|1] --csv PATH\n",
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

static void add_full_bodies(struct reb_simulation* sim) {
    struct reb_particle sun;
    memset(&sun, 0, sizeof(sun));
    sun.m = 1.0;
    add_cartesian(sim, &sun);
    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        const struct reb_particle p = particle_from_state(&GALILEAN_BODIES[i]);
        add_cartesian(sim, &p);
    }
    reb_simulation_move_to_com(sim);
}

static int set_global_split_integrator(
    struct reb_simulation* sim,
    const char* name,
    int whfast_corrector,
    int whfast_corrector2
) {
    if (strcmp(name, "whfast") == 0) {
        if (reb_bridge_set_integrator_whfast(sim) != 0) return -1;
        sim->ri_whfast.corrector = (unsigned int)whfast_corrector;
        sim->ri_whfast.corrector2 = (unsigned int)whfast_corrector2;
        sim->ri_whfast.safe_mode = 0;
        sim->ri_whfast.keep_unsynchronized = 1;
        return 0;
    }
    if (strcmp(name, "saba") == 0) {
        sim->integrator = REB_INTEGRATOR_SABA;
        sim->ri_whfast.coordinates = REB_WHFAST_COORDINATES_JACOBI;
        sim->ri_saba.type = REB_SABA_10_6_4;
        sim->ri_saba.safe_mode = 1;
        return 0;
    }
    return -1;
}

static struct reb_simulation* make_global_split(
    double dt_days,
    const char* integrator_name,
    int whfast_corrector,
    int whfast_corrector2
) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (set_global_split_integrator(sim, integrator_name, whfast_corrector, whfast_corrector2) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    sim->dt = dt_days / 365.25;
    add_full_bodies(sim);
    return sim;
}

static struct reb_simulation* make_ias15(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_ias15(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    add_full_bodies(sim);
    return sim;
}

static void fill_bodies(const struct reb_simulation* sim, struct reb_particle bodies[BODY_COUNT]) {
    for (int i = 0; i < BODY_COUNT; i++) bodies[i] = sim->particles[i];
}

static void fill_split_bodies(struct reb_simulation* sim, struct reb_particle bodies[BODY_COUNT]) {
    reb_simulation_synchronize(sim);
    fill_bodies(sim, bodies);
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
    options->dt_days = 0.05;
    options->whfast_corrector = 0;
    options->whfast_corrector2 = 0;
    options->integrator_name = "whfast";
    options->csv_path = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--integrator") == 0) {
            if (++i >= argc) return -1;
            if (strcmp(argv[i], "whfast") != 0 && strcmp(argv[i], "saba") != 0) return -1;
            options->integrator_name = argv[i];
        } else if (strcmp(argv[i], "--years") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->years) != 0 || options->years <= 0.0) return -1;
        } else if (strcmp(argv[i], "--samples") == 0) {
            if (++i >= argc || parse_int_arg(argv[i], &options->samples) != 0) return -1;
        } else if (strcmp(argv[i], "--dt-days") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->dt_days) != 0 || options->dt_days <= 0.0) return -1;
        } else if (strcmp(argv[i], "--whfast-corrector") == 0) {
            if (++i >= argc || parse_nonnegative_int_arg(argv[i], &options->whfast_corrector) != 0) return -1;
            if (!valid_whfast_corrector(options->whfast_corrector)) return -1;
        } else if (strcmp(argv[i], "--whfast-corrector2") == 0) {
            if (++i >= argc || parse_nonnegative_int_arg(argv[i], &options->whfast_corrector2) != 0) return -1;
            if (options->whfast_corrector2 != 0 && options->whfast_corrector2 != 1) return -1;
        } else if (strcmp(argv[i], "--csv") == 0) {
            if (++i >= argc) return -1;
            options->csv_path = argv[i];
        } else {
            return -1;
        }
    }

    if (strcmp(options->integrator_name, "whfast") != 0 &&
        (options->whfast_corrector != 0 || options->whfast_corrector2 != 0)) {
        return -1;
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

    struct reb_simulation* split = make_global_split(
        options.dt_days,
        options.integrator_name,
        options.whfast_corrector,
        options.whfast_corrector2
    );
    struct reb_simulation* ias15 = make_ias15();
    if (!split || !ias15) {
        if (split) reb_simulation_free(split);
        if (ias15) reb_simulation_free(ias15);
        fclose(csv);
        fprintf(stderr, "failed to initialize simulations\n");
        return 1;
    }

    fprintf(
        csv,
        "time_yr,split_step_seconds,ias15_step_seconds,split_cumulative_seconds,ias15_cumulative_seconds,cumulative_speedup,"
        "whfast_corrector,whfast_corrector2,"
        "split_ok,ias15_ok,split_escaped,ias15_escaped,"
        "io_error_au,europa_error_au,ganymede_error_au,callisto_error_au,max_moon_error_au,"
        "split_laplace_deg,ias15_laplace_deg,laplace_diff_deg,"
        "split_period_ratio_io_europa,ias15_period_ratio_io_europa,period_ratio_io_europa_diff,"
        "split_period_ratio_europa_ganymede,ias15_period_ratio_europa_ganymede,period_ratio_europa_ganymede_diff,"
        "split_max_e,ias15_max_e,split_max_distance_au,ias15_max_distance_au\n"
    );

    printf("Galilean naive global-%s vs IAS15 pilot\n", options.integrator_name);
    printf("epoch=%s source=%s\n", GALILEAN_HORIZONS_EPOCH, GALILEAN_HORIZONS_SOURCE);
    printf("years=%.17g samples=%d dt_days=%.17g whfast_corrector=%d whfast_corrector2=%d\n",
        options.years,
        options.samples,
        options.dt_days,
        options.whfast_corrector,
        options.whfast_corrector2
    );

    double split_total_seconds = 0.0;
    double ias15_total_seconds = 0.0;
    double max_moon_error = 0.0;
    double max_laplace_diff_deg = 0.0;
    double max_period_ratio_io_europa_diff = 0.0;
    double max_period_ratio_europa_ganymede_diff = 0.0;
    int split_escaped = 0;
    int ias15_escaped = 0;

    for (int sample = 1; sample <= options.samples; sample++) {
        const double target = options.years * (double)sample / (double)options.samples;

        const double split_start = monotonic_seconds();
        if (reb_simulation_integrate(split, target) != REB_STATUS_SUCCESS) {
            fprintf(stderr, "%s integration failed at %.17g yr\n", options.integrator_name, target);
            reb_simulation_free(split);
            reb_simulation_free(ias15);
            fclose(csv);
            return 1;
        }
        const double split_step_seconds = monotonic_seconds() - split_start;
        split_total_seconds += split_step_seconds;

        const double ias15_start = monotonic_seconds();
        if (reb_simulation_integrate(ias15, target) != REB_STATUS_SUCCESS) {
            fprintf(stderr, "ias15 integration failed at %.17g yr\n", target);
            reb_simulation_free(split);
            reb_simulation_free(ias15);
            fclose(csv);
            return 1;
        }
        const double ias15_step_seconds = monotonic_seconds() - ias15_start;
        ias15_total_seconds += ias15_step_seconds;

        struct reb_particle split_bodies[BODY_COUNT];
        struct reb_particle ias15_bodies[BODY_COUNT];
        fill_split_bodies(split, split_bodies);
        fill_bodies(ias15, ias15_bodies);

        const struct system_metrics ias15_metrics = compute_metrics(ias15_bodies, NULL, ias15->G);
        const struct system_metrics split_metrics = compute_metrics(split_bodies, ias15_bodies, split->G);
        const double laplace_diff_deg = fabs(angle_diff_rad(
            split_metrics.laplace_deg * M_PI / 180.0,
            ias15_metrics.laplace_deg * M_PI / 180.0
        )) * 180.0 / M_PI;
        const double period_ratio_io_europa_diff =
            fabs(split_metrics.period_ratio_io_europa - ias15_metrics.period_ratio_io_europa);
        const double period_ratio_europa_ganymede_diff =
            fabs(split_metrics.period_ratio_europa_ganymede - ias15_metrics.period_ratio_europa_ganymede);

        if (split_metrics.max_moon_error > max_moon_error) max_moon_error = split_metrics.max_moon_error;
        if (laplace_diff_deg > max_laplace_diff_deg) max_laplace_diff_deg = laplace_diff_deg;
        if (period_ratio_io_europa_diff > max_period_ratio_io_europa_diff) {
            max_period_ratio_io_europa_diff = period_ratio_io_europa_diff;
        }
        if (period_ratio_europa_ganymede_diff > max_period_ratio_europa_ganymede_diff) {
            max_period_ratio_europa_ganymede_diff = period_ratio_europa_ganymede_diff;
        }
        if (split_metrics.escaped) split_escaped = 1;
        if (ias15_metrics.escaped) ias15_escaped = 1;

        fprintf(
            csv,
            "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%d,%d,"
            "%d,%d,%d,%d,"
            "%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,%.17g\n",
            target,
            split_step_seconds,
            ias15_step_seconds,
            split_total_seconds,
            ias15_total_seconds,
            split_total_seconds > 0.0 ? ias15_total_seconds / split_total_seconds : 0.0,
            options.whfast_corrector,
            options.whfast_corrector2,
            split_metrics.ok,
            ias15_metrics.ok,
            split_metrics.escaped,
            ias15_metrics.escaped,
            split_metrics.moon_error[0],
            split_metrics.moon_error[1],
            split_metrics.moon_error[2],
            split_metrics.moon_error[3],
            split_metrics.max_moon_error,
            split_metrics.laplace_deg,
            ias15_metrics.laplace_deg,
            laplace_diff_deg,
            split_metrics.period_ratio_io_europa,
            ias15_metrics.period_ratio_io_europa,
            period_ratio_io_europa_diff,
            split_metrics.period_ratio_europa_ganymede,
            ias15_metrics.period_ratio_europa_ganymede,
            period_ratio_europa_ganymede_diff,
            split_metrics.max_eccentricity,
            ias15_metrics.max_eccentricity,
            split_metrics.max_distance_au,
            ias15_metrics.max_distance_au
        );
        fflush(csv);
    }

    printf("%s_total_seconds=%.9f\n", options.integrator_name, split_total_seconds);
    printf("ias15_total_seconds=%.9f\n", ias15_total_seconds);
    printf("ias15_over_%s_speedup=%.9f\n", options.integrator_name, split_total_seconds > 0.0 ? ias15_total_seconds / split_total_seconds : 0.0);
    printf("max_moon_error_au=%.17g\n", max_moon_error);
    printf("max_laplace_diff_deg=%.17g\n", max_laplace_diff_deg);
    printf("max_period_ratio_io_europa_diff=%.17g\n", max_period_ratio_io_europa_diff);
    printf("max_period_ratio_europa_ganymede_diff=%.17g\n", max_period_ratio_europa_ganymede_diff);
    printf("%s_escaped=%d\n", options.integrator_name, split_escaped);
    printf("ias15_escaped=%d\n", ias15_escaped);
    printf("wrote %s\n", options.csv_path);

    reb_simulation_free(split);
    reb_simulation_free(ias15);
    fclose(csv);
    return 0;
}
