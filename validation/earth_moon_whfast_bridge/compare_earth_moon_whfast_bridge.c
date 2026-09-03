/* Earth-Moon and Galilean satellite Bridge/WHFast/IAS15 triangle comparison. */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

struct barycenters {
    struct reb_particle earth_moon;
    struct reb_particle jovian;
};

struct comparison_maxima {
    double moon_relative_error_au;
    double moon_phase_arc_error_au;
    double earth_barycenter_error_au;
    double jupiter_barycenter_error_au;
    double io_relative_error_au;
    double europa_relative_error_au;
    double ganymede_relative_error_au;
};

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

static double pair_distance(const struct reb_particle* primary, const struct reb_particle* secondary) {
    return distance_particles(primary, secondary);
}

static double pair_phase_rad(const struct reb_particle* primary, const struct reb_particle* secondary) {
    return atan2(secondary->y - primary->y, secondary->x - primary->x);
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
    sim->exact_finish_time = 1;
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

static struct reb_simulation* make_sub_earth_moon(void) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    add_cartesian(sim, 3.0034896e-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    add_cartesian(sim, 3.6940e-8, 9.53035527126304600e-04, 2.17512105974572998e-03, 2.09676507944608496e-04, -2.09285374018137460e-01, 8.77025716574588027e-02, 1.16388188872971468e-03);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_sub_jovian(void) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    reb_simulation_add_fmt(sim, "m", 9.547919e-4);
    add_cartesian(sim, 2.388e-8, 2.48120485765584205e-03, -1.34556438159649208e-03, -1.31630007347439604e-05, 1.75341752994026812e+00, 3.20163492363371871e+00, 1.39542430489110164e-01);
    add_cartesian(sim, 2.527e-8, 5.51813635549570501e-04, -4.47551321186325738e-03, -1.27610904592036212e-04, 2.85755306406471643e+00, 3.72439324002732552e-01, 7.16871253077594472e-02);
    add_cartesian(sim, 7.804e-8, 6.76397567656043287e-03, -2.29852901298470890e-03, 1.01873701583559195e-05, 7.40898140660848692e-01, 2.17447712709783314e+00, 9.37799108604124437e-02);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct bridge_state make_bridge_state(double dt_outer_days, int earth_moon_ratio, int jovian_ratio, enum reb_bridge_method method) {
    struct bridge_state state;
    memset(&state, 0, sizeof(state));
    const double dt_outer = dt_outer_days / 365.25;
    state.main_sim = make_main_solar_system();
    state.earth_moon = make_sub_earth_moon();
    state.jovian = make_sub_jovian();
    if (!state.main_sim || !state.earth_moon || !state.jovian) return state;
    state.bridge = reb_bridge_create(state.main_sim, dt_outer);
    if (!state.bridge) return state;
    reb_bridge_set_owns_main(state.bridge, 1);
    if (reb_bridge_set_method(state.bridge, method) != 0) return state;
    if (reb_bridge_add_subsystem(state.bridge, MAIN_EARTH_MOON_HOST, state.earth_moon, dt_outer / (double)earth_moon_ratio, 1) != 0) return state;
    if (reb_bridge_add_subsystem(state.bridge, MAIN_JUPITER_HOST, state.jovian, dt_outer / (double)jovian_ratio, 1) != 0) return state;
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

static void fill_bridge_bodies(const struct bridge_state* state, struct reb_particle bodies[BODY_COUNT]) {
    const struct reb_particle* host_em = &state->main_sim->particles[MAIN_EARTH_MOON_HOST];
    const struct reb_particle* host_jov = &state->main_sim->particles[MAIN_JUPITER_HOST];
    bodies[BODY_SUN] = state->main_sim->particles[MAIN_SUN];
    bodies[BODY_MERCURY] = state->main_sim->particles[MAIN_MERCURY];
    bodies[BODY_VENUS] = state->main_sim->particles[MAIN_VENUS];
    bodies[BODY_EARTH] = val_particle_add(host_em, &state->earth_moon->particles[EM_EARTH]);
    bodies[BODY_MOON] = val_particle_add(host_em, &state->earth_moon->particles[EM_MOON]);
    bodies[BODY_MARS] = state->main_sim->particles[MAIN_MARS];
    bodies[BODY_JUPITER] = val_particle_add(host_jov, &state->jovian->particles[JOVIAN_JUPITER]);
    bodies[BODY_IO] = val_particle_add(host_jov, &state->jovian->particles[JOVIAN_IO]);
    bodies[BODY_EUROPA] = val_particle_add(host_jov, &state->jovian->particles[JOVIAN_EUROPA]);
    bodies[BODY_GANYMEDE] = val_particle_add(host_jov, &state->jovian->particles[JOVIAN_GANYMEDE]);
    bodies[BODY_SATURN] = state->main_sim->particles[MAIN_SATURN];
    bodies[BODY_URANUS] = state->main_sim->particles[MAIN_URANUS];
    bodies[BODY_NEPTUNE] = state->main_sim->particles[MAIN_NEPTUNE];
}

static void fill_whfast_bodies(const struct reb_simulation* sim, struct reb_particle bodies[BODY_COUNT]) {
    for (int i = 0; i < BODY_COUNT; i++) {
        bodies[i] = sim->particles[i];
    }
}

static struct barycenters compute_barycenters(const struct reb_particle bodies[BODY_COUNT]) {
    struct barycenters out;
    const int em_ids[] = { BODY_EARTH, BODY_MOON };
    const int jov_ids[] = { BODY_JUPITER, BODY_IO, BODY_EUROPA, BODY_GANYMEDE };
    out.earth_moon = val_barycenter(bodies, em_ids, 2);
    out.jovian = val_barycenter(bodies, jov_ids, 4);
    return out;
}

static void update_comparison_maxima(
    struct comparison_maxima* maxima,
    const struct reb_particle a_bodies[BODY_COUNT],
    const struct reb_particle b_bodies[BODY_COUNT]
) {
    const struct barycenters a_bc = compute_barycenters(a_bodies);
    const struct barycenters b_bc = compute_barycenters(b_bodies);

    const double moon_relative_error_au = relative_position_delta(
        &a_bodies[BODY_EARTH],
        &a_bodies[BODY_MOON],
        &b_bodies[BODY_EARTH],
        &b_bodies[BODY_MOON]
    );
    const double moon_phase_arc_error_au =
        fabs(angle_diff_rad(
            pair_phase_rad(&a_bodies[BODY_EARTH], &a_bodies[BODY_MOON]),
            pair_phase_rad(&b_bodies[BODY_EARTH], &b_bodies[BODY_MOON])
        )) *
        0.5 * (
            pair_distance(&a_bodies[BODY_EARTH], &a_bodies[BODY_MOON]) +
            pair_distance(&b_bodies[BODY_EARTH], &b_bodies[BODY_MOON])
        );
    const double earth_barycenter_error_au = distance_particles(&a_bc.earth_moon, &b_bc.earth_moon);
    const double jupiter_barycenter_error_au = distance_particles(&a_bc.jovian, &b_bc.jovian);
    const double io_relative_error_au = relative_position_delta(&a_bodies[BODY_JUPITER], &a_bodies[BODY_IO], &b_bodies[BODY_JUPITER], &b_bodies[BODY_IO]);
    const double europa_relative_error_au = relative_position_delta(&a_bodies[BODY_JUPITER], &a_bodies[BODY_EUROPA], &b_bodies[BODY_JUPITER], &b_bodies[BODY_EUROPA]);
    const double ganymede_relative_error_au = relative_position_delta(&a_bodies[BODY_JUPITER], &a_bodies[BODY_GANYMEDE], &b_bodies[BODY_JUPITER], &b_bodies[BODY_GANYMEDE]);

    if (moon_relative_error_au > maxima->moon_relative_error_au) maxima->moon_relative_error_au = moon_relative_error_au;
    if (moon_phase_arc_error_au > maxima->moon_phase_arc_error_au) maxima->moon_phase_arc_error_au = moon_phase_arc_error_au;
    if (earth_barycenter_error_au > maxima->earth_barycenter_error_au) maxima->earth_barycenter_error_au = earth_barycenter_error_au;
    if (jupiter_barycenter_error_au > maxima->jupiter_barycenter_error_au) maxima->jupiter_barycenter_error_au = jupiter_barycenter_error_au;
    if (io_relative_error_au > maxima->io_relative_error_au) maxima->io_relative_error_au = io_relative_error_au;
    if (europa_relative_error_au > maxima->europa_relative_error_au) maxima->europa_relative_error_au = europa_relative_error_au;
    if (ganymede_relative_error_au > maxima->ganymede_relative_error_au) maxima->ganymede_relative_error_au = ganymede_relative_error_au;
}


static void print_usage(const char* argv0) {
    fprintf(
        stderr,
        "usage: %s [--years N] [--samples N] [--dt-outer-days N] [--earth-moon-ratio N] [--jovian-ratio N] [--bridge-method order2|eos-lf4-2|eos-lf8-6-4] [--csv PATH]\n",
        argv0
    );
}

int main(int argc, char** argv) {
    double years = 100.0;
    int samples = 101;
    double dt_outer_days = 1.0;
    int earth_moon_ratio = 20;
    int jovian_ratio = 80;
    enum reb_bridge_method method = REB_BRIDGE_METHOD_EOS_LF4_2;
    const char* csv_path = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--years") == 0) {
            if (i + 1 >= argc || val_parse_double(argv[++i], &years) != 0 || years <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--samples") == 0) {
            if (i + 1 >= argc || val_parse_int(argv[++i], &samples) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--dt-outer-days") == 0) {
            if (i + 1 >= argc || val_parse_double(argv[++i], &dt_outer_days) != 0 || dt_outer_days <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--earth-moon-ratio") == 0) {
            if (i + 1 >= argc || val_parse_int(argv[++i], &earth_moon_ratio) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--jovian-ratio") == 0) {
            if (i + 1 >= argc || val_parse_int(argv[++i], &jovian_ratio) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--bridge-method") == 0) {
            if (i + 1 >= argc || val_parse_bridge_method(argv[++i], &method) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--csv") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            csv_path = argv[++i];
        } else {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    const double whfast_dt_days = dt_outer_days / (double)(earth_moon_ratio > jovian_ratio ? earth_moon_ratio : jovian_ratio);
    struct bridge_state bridge_state = make_bridge_state(dt_outer_days, earth_moon_ratio, jovian_ratio, method);
    struct reb_simulation* whfast = ss_make_whfast(whfast_dt_days);
    struct reb_simulation* ias15 = ss_make_ias15();
    if (!bridge_state.bridge || !bridge_state.main_sim || !bridge_state.earth_moon || !bridge_state.jovian || !whfast || !ias15) {
        free_bridge_state(&bridge_state);
        if (whfast) reb_simulation_free(whfast);
        if (ias15) reb_simulation_free(ias15);
        return EXIT_FAILURE;
    }

    double bridge_total_seconds = 0.0;
    double whfast_total_seconds = 0.0;
    double ias15_total_seconds = 0.0;
    struct comparison_maxima bridge_vs_ias15;
    struct comparison_maxima whfast_vs_ias15;
    struct comparison_maxima bridge_vs_whfast;
    memset(&bridge_vs_ias15, 0, sizeof(bridge_vs_ias15));
    memset(&whfast_vs_ias15, 0, sizeof(whfast_vs_ias15));
    memset(&bridge_vs_whfast, 0, sizeof(bridge_vs_whfast));

    for (int sample = 1; sample <= samples; sample++) {
        const double target = years * (double)sample / (double)samples;
        double bridge_duration = target - bridge_state.main_sim->t;
        const double bridge_start = monotonic_seconds();
        /* Advance in outer-sized chunks; finish the final partial chunk with
         * reb_bridge_step so 100 yr need not be an integer number of 7-day steps. */
        const double outer_step = dt_outer_days / 365.25;
        while (bridge_duration > outer_step * 1.0000000001) {
            if (reb_bridge_step(bridge_state.bridge, outer_step) != 0) break;
            bridge_duration -= outer_step;
        }
        if (bridge_duration > 0.0 && reb_bridge_step(bridge_state.bridge, bridge_duration) != 0) {
            free_bridge_state(&bridge_state);
            reb_simulation_free(whfast);
            reb_simulation_free(ias15);
            return EXIT_FAILURE;
        }
        bridge_total_seconds += monotonic_seconds() - bridge_start;

        const double whfast_start = monotonic_seconds();
        if (reb_simulation_integrate(whfast, target) != REB_STATUS_SUCCESS) {
            free_bridge_state(&bridge_state);
            reb_simulation_free(whfast);
            reb_simulation_free(ias15);
            return EXIT_FAILURE;
        }
        whfast_total_seconds += monotonic_seconds() - whfast_start;

        const double ias15_start = monotonic_seconds();
        if (reb_simulation_integrate(ias15, target) != REB_STATUS_SUCCESS) {
            free_bridge_state(&bridge_state);
            reb_simulation_free(whfast);
            reb_simulation_free(ias15);
            return EXIT_FAILURE;
        }
        ias15_total_seconds += monotonic_seconds() - ias15_start;

        struct reb_particle bridge_bodies[BODY_COUNT];
        struct reb_particle whfast_bodies[BODY_COUNT];
        struct reb_particle ias15_bodies[BODY_COUNT];
        fill_bridge_bodies(&bridge_state, bridge_bodies);
        fill_whfast_bodies(whfast, whfast_bodies);
        fill_whfast_bodies(ias15, ias15_bodies);

        update_comparison_maxima(&bridge_vs_ias15, bridge_bodies, ias15_bodies);
        update_comparison_maxima(&whfast_vs_ias15, whfast_bodies, ias15_bodies);
        update_comparison_maxima(&bridge_vs_whfast, bridge_bodies, whfast_bodies);
    }

    printf("years=%.9f\n", years);
    printf("samples=%d\n", samples);
    printf("dt_outer_days=%.9f\n", dt_outer_days);
    printf("earth_moon_ratio=%d\n", earth_moon_ratio);
    printf("jovian_ratio=%d\n", jovian_ratio);
    printf("bridge_method=%s\n", val_bridge_method_name(method));
    printf("whfast_dt_days=%.9f\n", whfast_dt_days);
    printf("bridge_total_seconds=%.9f\n", bridge_total_seconds);
    printf("whfast_total_seconds=%.9f\n", whfast_total_seconds);
    printf("ias15_total_seconds=%.9f\n", ias15_total_seconds);
    printf("whfast_over_bridge_speedup=%.9f\n", bridge_total_seconds > 0.0 ? whfast_total_seconds / bridge_total_seconds : 0.0);
    printf("ias15_over_bridge_speedup=%.9f\n", bridge_total_seconds > 0.0 ? ias15_total_seconds / bridge_total_seconds : 0.0);
    printf("ias15_over_whfast_speedup=%.9f\n", whfast_total_seconds > 0.0 ? ias15_total_seconds / whfast_total_seconds : 0.0);
    printf("bridge_vs_ias15_moon_relative_error_au=%.17g\n", bridge_vs_ias15.moon_relative_error_au);
    printf("whfast_vs_ias15_moon_relative_error_au=%.17g\n", whfast_vs_ias15.moon_relative_error_au);
    printf("bridge_vs_whfast_moon_relative_error_au=%.17g\n", bridge_vs_whfast.moon_relative_error_au);
    printf("bridge_vs_ias15_moon_phase_arc_error_au=%.17g\n", bridge_vs_ias15.moon_phase_arc_error_au);
    printf("whfast_vs_ias15_moon_phase_arc_error_au=%.17g\n", whfast_vs_ias15.moon_phase_arc_error_au);
    printf("bridge_vs_whfast_moon_phase_arc_error_au=%.17g\n", bridge_vs_whfast.moon_phase_arc_error_au);
    printf("bridge_vs_ias15_earth_barycenter_error_au=%.17g\n", bridge_vs_ias15.earth_barycenter_error_au);
    printf("whfast_vs_ias15_earth_barycenter_error_au=%.17g\n", whfast_vs_ias15.earth_barycenter_error_au);
    printf("bridge_vs_whfast_earth_barycenter_error_au=%.17g\n", bridge_vs_whfast.earth_barycenter_error_au);
    printf("bridge_vs_ias15_jupiter_barycenter_error_au=%.17g\n", bridge_vs_ias15.jupiter_barycenter_error_au);
    printf("whfast_vs_ias15_jupiter_barycenter_error_au=%.17g\n", whfast_vs_ias15.jupiter_barycenter_error_au);
    printf("bridge_vs_whfast_jupiter_barycenter_error_au=%.17g\n", bridge_vs_whfast.jupiter_barycenter_error_au);

    if (csv_path) {
        FILE* fp = fopen(csv_path, "a");
        if (!fp) {
            free_bridge_state(&bridge_state);
            reb_simulation_free(whfast);
            reb_simulation_free(ias15);
            return EXIT_FAILURE;
        }
        fprintf(
            fp,
            "%s,%.17g,%d,%.17g,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            val_bridge_method_name(method),
            years,
            samples,
            dt_outer_days,
            earth_moon_ratio,
            jovian_ratio,
            whfast_dt_days,
            bridge_total_seconds,
            whfast_total_seconds,
            ias15_total_seconds,
            bridge_total_seconds > 0.0 ? whfast_total_seconds / bridge_total_seconds : 0.0,
            bridge_total_seconds > 0.0 ? ias15_total_seconds / bridge_total_seconds : 0.0,
            whfast_total_seconds > 0.0 ? ias15_total_seconds / whfast_total_seconds : 0.0,
            bridge_vs_ias15.moon_relative_error_au,
            whfast_vs_ias15.moon_relative_error_au,
            bridge_vs_whfast.moon_relative_error_au,
            bridge_vs_ias15.moon_phase_arc_error_au,
            whfast_vs_ias15.moon_phase_arc_error_au,
            bridge_vs_whfast.moon_phase_arc_error_au,
            bridge_vs_ias15.earth_barycenter_error_au,
            whfast_vs_ias15.earth_barycenter_error_au,
            bridge_vs_whfast.earth_barycenter_error_au,
            bridge_vs_ias15.jupiter_barycenter_error_au,
            whfast_vs_ias15.jupiter_barycenter_error_au,
            bridge_vs_whfast.jupiter_barycenter_error_au,
            bridge_vs_ias15.io_relative_error_au,
            whfast_vs_ias15.io_relative_error_au,
            bridge_vs_whfast.io_relative_error_au,
            bridge_vs_ias15.europa_relative_error_au,
            whfast_vs_ias15.europa_relative_error_au,
            bridge_vs_whfast.europa_relative_error_au,
            bridge_vs_ias15.ganymede_relative_error_au,
            whfast_vs_ias15.ganymede_relative_error_au,
            bridge_vs_whfast.ganymede_relative_error_au
        );
        fclose(fp);
    }

    free_bridge_state(&bridge_state);
    reb_simulation_free(whfast);
    reb_simulation_free(ias15);
    return EXIT_SUCCESS;
}
