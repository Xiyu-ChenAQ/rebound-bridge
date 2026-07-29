/* Long-run Galilean energy benchmark: Bridge vs global WHFast, no IAS15. */

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
    BODY_SUN = 0,
    BODY_JUPITER = 1,
    BODY_IO = 2,
    BODY_EUROPA = 3,
    BODY_GANYMEDE = 4,
    BODY_CALLISTO = 5,
    BODY_COUNT = 6
};

struct run_options {
    const char* method_name;
    double years;
    double dt_outer_days;
    double dt_days;
    int inner_ratio;
    int realization_start;
    int realizations;
    double perturbation;
    int outputs_per_decade;
    int whfast_corrector;
    int whfast_corrector2;
    const char* csv_path;
};

struct bridge_state {
    struct reb_simulation* main_sim;
    struct reb_simulation* sub_sim;
    struct reb_bridge* bridge;
};

struct stability_metrics {
    int ok;
    int escaped;
    double max_eccentricity;
    double max_distance_au;
    double laplace_deg;
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

static void print_usage(const char* argv0) {
    fprintf(
        stderr,
        "usage: %s --method bridge|whfast [--years N] [--dt-outer-days N] [--dt-days N] "
        "[--inner-ratio N] [--realization-start N] [--realizations N] [--perturbation N] [--outputs-per-decade N] "
        "[--whfast-corrector N] [--whfast-corrector2 0|1] --csv PATH\n",
        argv0
    );
}

static uint64_t rng_next(uint64_t* state) {
    *state = (*state * 6364136223846793005ULL) + 1442695040888963407ULL;
    return *state;
}

static double uniform_unit(uint64_t* state) {
    const uint64_t x = rng_next(state) >> 11;
    return (double)x * (1.0 / 9007199254740992.0);
}

static double uniform_signed(uint64_t* state) {
    return 2.0 * uniform_unit(state) - 1.0;
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

static void perturb_phase_space(struct reb_particle* p, uint64_t* rng, double perturbation) {
    if (perturbation <= 0.0) return;
    p->x *= 1.0 + perturbation * uniform_signed(rng);
    p->y *= 1.0 + perturbation * uniform_signed(rng);
    p->z *= 1.0 + perturbation * uniform_signed(rng);
    p->vx *= 1.0 + perturbation * uniform_signed(rng);
    p->vy *= 1.0 + perturbation * uniform_signed(rng);
    p->vz *= 1.0 + perturbation * uniform_signed(rng);
}

static void add_cartesian(struct reb_simulation* sim, const struct reb_particle* p) {
    reb_simulation_add_fmt(sim, "m x y z vx vy vz", p->m, p->x, p->y, p->z, p->vx, p->vy, p->vz);
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

static void make_initial_bodies(int realization, double perturbation, struct reb_particle bodies[BODY_COUNT]) {
    struct reb_particle template_bodies[BODY_COUNT];
    memset(template_bodies, 0, sizeof(template_bodies));
    template_bodies[BODY_SUN].m = 1.0;
    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        template_bodies[i + 1] = particle_from_state(&GALILEAN_BODIES[i]);
    }

    const struct reb_particle template_host = host_com(template_bodies);
    struct reb_particle perturbed_host = template_host;
    uint64_t rng = 0x9e3779b97f4a7c15ULL ^ (uint64_t)(realization + 1) * 0xbf58476d1ce4e5b9ULL;

    memset(bodies, 0, sizeof(struct reb_particle) * BODY_COUNT);
    bodies[BODY_SUN] = template_bodies[BODY_SUN];
    if (perturbation > 0.0) {
        bodies[BODY_SUN].m *= 1.0 + perturbation * uniform_signed(&rng);
        perturb_phase_space(&perturbed_host, &rng, perturbation);
    }

    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        struct reb_particle relative = subtract_particle(&template_bodies[i + 1], &template_host);
        if (perturbation > 0.0) {
            relative.m *= 1.0 + perturbation * uniform_signed(&rng);
            perturb_phase_space(&relative, &rng, perturbation);
        }
        bodies[i + 1] = add_particle(&relative, &perturbed_host);
    }
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
            energy -= G * bodies[i].m * bodies[j].m / distance_particles(&bodies[i], &bodies[j]);
        }
    }
    return energy;
}

static double total_angular_momentum_norm(const struct reb_particle bodies[BODY_COUNT]) {
    double lx = 0.0;
    double ly = 0.0;
    double lz = 0.0;
    for (int i = 0; i < BODY_COUNT; i++) {
        const struct reb_particle* p = &bodies[i];
        lx += p->m * (p->y * p->vz - p->z * p->vy);
        ly += p->m * (p->z * p->vx - p->x * p->vz);
        lz += p->m * (p->x * p->vy - p->y * p->vx);
    }
    return norm3(lx, ly, lz);
}

static double total_linear_momentum_norm(const struct reb_particle bodies[BODY_COUNT]) {
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    for (int i = 0; i < BODY_COUNT; i++) {
        px += bodies[i].m * bodies[i].vx;
        py += bodies[i].m * bodies[i].vy;
        pz += bodies[i].m * bodies[i].vz;
    }
    return norm3(px, py, pz);
}

static double center_of_mass_distance(const struct reb_particle bodies[BODY_COUNT]) {
    double total_mass = 0.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    for (int i = 0; i < BODY_COUNT; i++) {
        total_mass += bodies[i].m;
        x += bodies[i].m * bodies[i].x;
        y += bodies[i].m * bodies[i].y;
        z += bodies[i].m * bodies[i].z;
    }
    return total_mass > 0.0 ? norm3(x, y, z) / total_mass : 0.0;
}

static struct stability_metrics compute_stability_metrics(const struct reb_particle bodies[BODY_COUNT], double G) {
    struct stability_metrics out;
    memset(&out, 0, sizeof(out));
    out.ok = 1;

    struct reb_orbit orbits[4];
    for (int i = 0; i < 4; i++) {
        int err = 0;
        const int body_id = BODY_IO + i;
        orbits[i] = reb_orbit_from_particle_err(G, bodies[body_id], bodies[BODY_JUPITER], &err);
        const double distance = distance_particles(&bodies[BODY_JUPITER], &bodies[body_id]);
        if (distance > out.max_distance_au) out.max_distance_au = distance;
        if (err || !isfinite(orbits[i].e) || !isfinite(orbits[i].l)) {
            out.ok = 0;
            out.escaped = 1;
            continue;
        }
        if (orbits[i].e > out.max_eccentricity) out.max_eccentricity = orbits[i].e;
        if (distance < R_JUPITER_AU || orbits[i].e >= 1.0) out.escaped = 1;
    }

    if (out.ok) {
        out.laplace_deg = wrap_rad(orbits[0].l - 3.0 * orbits[1].l + 2.0 * orbits[2].l) * 180.0 / M_PI;
    }
    return out;
}

static void fill_bridge_bodies(const struct bridge_state* state, struct reb_particle bodies[BODY_COUNT]) {
    const struct reb_particle host = state->main_sim->particles[1];
    bodies[BODY_SUN] = state->main_sim->particles[0];
    for (int i = 0; i < GALILEAN_BODY_COUNT; i++) {
        bodies[i + 1] = add_particle(&state->sub_sim->particles[i], &host);
    }
}

static void fill_sim_bodies(const struct reb_simulation* sim, struct reb_particle bodies[BODY_COUNT]) {
    for (int i = 0; i < BODY_COUNT; i++) bodies[i] = sim->particles[i];
}

static void fill_whfast_bodies(struct reb_simulation* sim, struct reb_particle bodies[BODY_COUNT]) {
    reb_simulation_synchronize(sim);
    fill_sim_bodies(sim, bodies);
}

static struct reb_simulation* make_whfast_sim(
    const struct reb_particle bodies[BODY_COUNT],
    double dt_days,
    int corrector,
    int corrector2
) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 0;
    sim->dt = dt_days / 365.25;
    sim->ri_whfast.corrector = (unsigned int)corrector;
    sim->ri_whfast.corrector2 = (unsigned int)corrector2;
    sim->ri_whfast.safe_mode = 0;
    sim->ri_whfast.keep_unsynchronized = 1;
    for (int i = 0; i < BODY_COUNT; i++) add_cartesian(sim, &bodies[i]);
    reb_simulation_move_to_com(sim);
    return sim;
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

    add_cartesian(state.main_sim, &bodies[BODY_SUN]);
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
    options->method_name = NULL;
    options->years = 1000.0;
    options->dt_outer_days = 0.5;
    options->dt_days = 0.5;
    options->inner_ratio = 40;
    options->realization_start = 0;
    options->realizations = 1;
    options->perturbation = 0.001;
    options->outputs_per_decade = 20;
    options->whfast_corrector = 0;
    options->whfast_corrector2 = 0;
    options->csv_path = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--method") == 0) {
            if (++i >= argc) return -1;
            if (strcmp(argv[i], "bridge") != 0 && strcmp(argv[i], "whfast") != 0) return -1;
            options->method_name = argv[i];
        } else if (strcmp(argv[i], "--years") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->years) != 0 || options->years <= 0.0) return -1;
        } else if (strcmp(argv[i], "--dt-outer-days") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->dt_outer_days) != 0 || options->dt_outer_days <= 0.0) return -1;
        } else if (strcmp(argv[i], "--dt-days") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->dt_days) != 0 || options->dt_days <= 0.0) return -1;
        } else if (strcmp(argv[i], "--inner-ratio") == 0) {
            if (++i >= argc || parse_int_arg(argv[i], &options->inner_ratio) != 0) return -1;
        } else if (strcmp(argv[i], "--realization-start") == 0) {
            if (++i >= argc || parse_nonnegative_int_arg(argv[i], &options->realization_start) != 0) return -1;
        } else if (strcmp(argv[i], "--realizations") == 0) {
            if (++i >= argc || parse_int_arg(argv[i], &options->realizations) != 0) return -1;
        } else if (strcmp(argv[i], "--perturbation") == 0) {
            if (++i >= argc || parse_double_arg(argv[i], &options->perturbation) != 0 || options->perturbation < 0.0) return -1;
        } else if (strcmp(argv[i], "--outputs-per-decade") == 0) {
            if (++i >= argc || parse_int_arg(argv[i], &options->outputs_per_decade) != 0) return -1;
        } else if (strcmp(argv[i], "--whfast-corrector") == 0) {
            if (++i >= argc || parse_nonnegative_int_arg(argv[i], &options->whfast_corrector) != 0) return -1;
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
    return options->method_name && options->csv_path ? 0 : -1;
}

static uint64_t next_log_output_step(uint64_t previous, uint64_t total_steps, int outputs_per_decade) {
    if (previous == 0) return 1;
    const double next = ceil(pow(10.0, log10((double)previous) + 1.0 / (double)outputs_per_decade));
    uint64_t out = (uint64_t)next;
    if (out <= previous) out = previous + 1;
    if (out > total_steps) out = total_steps;
    return out;
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
    fprintf(
        csv,
        "method,realization,step,time_yr,wall_seconds,energy,energy_rel,energy_rel_abs,"
        "angular_momentum,angular_momentum_rel,angular_momentum_rel_abs,linear_momentum_norm,com_distance_au,"
        "dt_outer_days,dt_days,inner_ratio,perturbation,whfast_corrector,whfast_corrector2,"
        "stability_ok,escaped,max_eccentricity,max_distance_au,laplace_deg\n"
    );

    printf("Galilean long energy benchmark\n");
    printf("method=%s years=%.17g dt_outer_days=%.17g dt_days=%.17g inner_ratio=%d realizations=%d perturbation=%.17g whfast_corrector=%d whfast_corrector2=%d\n",
        options.method_name,
        options.years,
        options.dt_outer_days,
        options.dt_days,
        options.inner_ratio,
        options.realizations,
        options.perturbation,
        options.whfast_corrector,
        options.whfast_corrector2
    );

    double total_wall_seconds = 0.0;
    double max_energy_rel_abs = 0.0;
    double max_angular_momentum_rel_abs = 0.0;

    for (int run = 0; run < options.realizations; run++) {
        const int realization = options.realization_start + run;
        struct reb_particle initial_bodies[BODY_COUNT];
        struct reb_particle current_bodies[BODY_COUNT];
        make_initial_bodies(realization, options.perturbation, initial_bodies);

        const double dt_days = strcmp(options.method_name, "bridge") == 0
            ? options.dt_outer_days
            : options.dt_days;
        const double dt_years = dt_days / 365.25;
        const uint64_t total_steps = (uint64_t)llround(options.years / dt_years);
        uint64_t next_output_step = 1;
        double initial_energy = 0.0;
        double initial_angular_momentum = 0.0;
        double start = monotonic_seconds();

        struct bridge_state bridge;
        memset(&bridge, 0, sizeof(bridge));
        struct reb_simulation* whfast = NULL;

        if (strcmp(options.method_name, "bridge") == 0) {
            bridge = make_bridge_state(initial_bodies, options.dt_outer_days, options.inner_ratio);
            if (!bridge.bridge || !bridge.main_sim || !bridge.sub_sim) {
                free_bridge_state(&bridge);
                fclose(csv);
                fprintf(stderr, "failed to initialize Bridge realization %d\n", realization);
                return 1;
            }
            fill_bridge_bodies(&bridge, current_bodies);
        } else {
            whfast = make_whfast_sim(
                initial_bodies,
                options.dt_days,
                options.whfast_corrector,
                options.whfast_corrector2
            );
            if (!whfast) {
                fclose(csv);
                fprintf(stderr, "failed to initialize WHFast realization %d\n", realization);
                return 1;
            }
            fill_whfast_bodies(whfast, current_bodies);
        }
        initial_energy = total_energy(current_bodies, REB_BRIDGE_G_AU_YR_MSUN);
        initial_angular_momentum = total_angular_momentum_norm(current_bodies);

        for (uint64_t step = 1; step <= total_steps; step++) {
            if (strcmp(options.method_name, "bridge") == 0) {
                if (reb_bridge_step(bridge.bridge, options.dt_outer_days / 365.25) != 0) {
                    free_bridge_state(&bridge);
                    fclose(csv);
                    fprintf(stderr, "Bridge step failed at realization %d step %llu\n", realization, (unsigned long long)step);
                    return 1;
                }
                if (step == next_output_step) fill_bridge_bodies(&bridge, current_bodies);
            } else {
                reb_simulation_step(whfast);
                if (step == next_output_step) fill_whfast_bodies(whfast, current_bodies);
            }

            if (step == next_output_step) {
                const double wall = monotonic_seconds() - start;
                const double energy = total_energy(current_bodies, REB_BRIDGE_G_AU_YR_MSUN);
                const double energy_rel = (energy - initial_energy) / fabs(initial_energy);
                const double energy_abs = fabs(energy_rel);
                const double angular_momentum = total_angular_momentum_norm(current_bodies);
                const double angular_momentum_rel =
                    (angular_momentum - initial_angular_momentum) / fabs(initial_angular_momentum);
                const double angular_momentum_abs = fabs(angular_momentum_rel);
                const double linear_momentum_norm = total_linear_momentum_norm(current_bodies);
                const double com_distance_au = center_of_mass_distance(current_bodies);
                const struct stability_metrics stability =
                    compute_stability_metrics(current_bodies, REB_BRIDGE_G_AU_YR_MSUN);
                if (energy_abs > max_energy_rel_abs) max_energy_rel_abs = energy_abs;
                if (angular_momentum_abs > max_angular_momentum_rel_abs) {
                    max_angular_momentum_rel_abs = angular_momentum_abs;
                }
                fprintf(
                    csv,
                    "%s,%d,%llu,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%d,%.17g,%d,%d,%d,%d,%.17g,%.17g,%.17g\n",
                    options.method_name,
                    realization,
                    (unsigned long long)step,
                    (double)step * dt_years,
                    wall,
                    energy,
                    energy_rel,
                    energy_abs,
                    angular_momentum,
                    angular_momentum_rel,
                    angular_momentum_abs,
                    linear_momentum_norm,
                    com_distance_au,
                    options.dt_outer_days,
                    options.dt_days,
                    options.inner_ratio,
                    options.perturbation,
                    options.whfast_corrector,
                    options.whfast_corrector2,
                    stability.ok,
                    stability.escaped,
                    stability.max_eccentricity,
                    stability.max_distance_au,
                    stability.laplace_deg
                );
                fflush(csv);
                next_output_step = next_log_output_step(step, total_steps, options.outputs_per_decade);
            }
        }

        const double elapsed = monotonic_seconds() - start;
        total_wall_seconds += elapsed;
        if (bridge.bridge) free_bridge_state(&bridge);
        if (whfast) reb_simulation_free(whfast);
        printf("realization=%d wall_seconds=%.9f\n", realization, elapsed);
    }

    printf("total_wall_seconds=%.9f\n", total_wall_seconds);
    printf("max_energy_rel_abs=%.17g\n", max_energy_rel_abs);
    printf("max_angular_momentum_rel_abs=%.17g\n", max_angular_momentum_rel_abs);
    printf("wrote %s\n", options.csv_path);
    fclose(csv);
    return 0;
}
