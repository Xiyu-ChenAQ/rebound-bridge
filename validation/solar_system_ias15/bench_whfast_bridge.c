/* SPDX-License-Identifier: GPL-3.0-only */

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

#include "rebound_bridge.h"

enum main_body {
    MAIN_EARTH_MOON_HOST = 3,
    MAIN_JUPITER_HOST = 5
};

struct bridge_state {
    struct reb_simulation* main_sim;
    struct reb_simulation* earth_moon;
    struct reb_simulation* jovian;
    struct reb_bridge* bridge;
};

static double monotonic_seconds(void) {
#ifdef _WIN32
    static LARGE_INTEGER frequency;
    static int initialized = 0;
    LARGE_INTEGER counter;
    if (!initialized) {
        QueryPerformanceFrequency(&frequency);
        initialized = 1;
    }
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / (double)frequency.QuadPart;
#else
    struct timespec ts;
#if defined(CLOCK_MONOTONIC)
    clock_gettime(CLOCK_MONOTONIC, &ts);
#else
    timespec_get(&ts, TIME_UTC);
#endif
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
#endif
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

static void add_full_physical_bodies(struct reb_simulation* sim) {
    reb_simulation_add_fmt(sim, "m", 1.0);
    add_cartesian(sim, 1.660120e-7, -2.15201304377667713e-01, -4.09207615700272087e-01, -1.37032695919236600e-02, 7.02443404629697010e+00, -4.28716799788450942e+00, -9.94632701667513053e-01);
    add_cartesian(sim, 2.4478383e-6, 8.88772496163565545e-02, -7.21762382317080675e-01, -1.50440565750619305e-02, 7.28255910865529810e+00, 8.76733642322772511e-01, -4.08156573658943866e-01);
    add_cartesian(sim, 3.0034896e-6, -1.74269758548305093e-01, 9.67785674324191270e-01, -5.68660813624550973e-05, -6.28653514426191418e+00, -1.13726197010495311e+00, 1.00945156455726687e-04);
    add_cartesian(sim, 3.6940e-8, -1.74269758548305093e-01 + 9.53035527126304600e-04, 9.67785674324191270e-01 + 2.17512105974572998e-03, -5.68660813624550973e-05 + 2.09676507944608496e-04, -6.28653514426191418e+00 - 2.09285374018137460e-01, -1.13726197010495311e+00 + 8.77025716574588027e-02, 1.00945156455726687e-04 + 1.16388188872971468e-03);
    add_cartesian(sim, 3.227151e-7, 3.40579676862215075e-01, -1.38700201594525407e+00, -3.74172267877010828e-02, 5.15711521620981550e+00, 1.65830921706158052e+00, -9.17089273335817273e-02);
    add_cartesian(sim, 9.547919e-4, -1.69400303018584597e+00, 4.92888271463623262e+00, 1.74262399033241998e-02, -2.64098245722763814e+00, -7.67587434003239921e-01, 6.22754635347801122e-02);
    add_cartesian(sim, 2.388e-8, -1.69400303018584597e+00 + 2.48120485765584205e-03, 4.92888271463623262e+00 - 1.34556438159649208e-03, 1.74262399033241998e-02 - 1.31630007347439604e-05, -2.64098245722763814e+00 + 1.75341752994026812e+00, -7.67587434003239921e-01 + 3.20163492363371871e+00, 6.22754635347801122e-02 + 1.39542430489110164e-01);
    add_cartesian(sim, 2.527e-8, -1.69400303018584597e+00 + 5.51813635549570501e-04, 4.92888271463623262e+00 - 4.47551321186325738e-03, 1.74262399033241998e-02 - 1.27610904592036212e-04, -2.64098245722763814e+00 + 2.85755306406471643e+00, -7.67587434003239921e-01 + 3.72439324002732552e-01, 6.22754635347801122e-02 + 7.16871253077594472e-02);
    add_cartesian(sim, 7.804e-8, -1.69400303018584597e+00 + 6.76397567656043287e-03, 4.92888271463623262e+00 - 2.29852901298470890e-03, 1.74262399033241998e-02 + 1.01873701583559195e-05, -2.64098245722763814e+00 + 7.40898140660848692e-01, -7.67587434003239921e-01 + 2.17447712709783314e+00, 6.22754635347801122e-02 + 9.37799108604124437e-02);
    add_cartesian(sim, 2.8588567e-4, 9.50734336522187107e+00, 2.57738457802152210e-01, -3.82935528507753087e-01, -1.69141631745304527e-01, 2.03249650346815258e+00, -2.86682018818476998e-02);
    add_cartesian(sim, 4.366244e-5, 9.88031797519427535e+00, 1.68000154098179486e+01, -6.57197991373143220e-02, -1.25149298434788658e+00, 6.61240423623731699e-01, 1.86840216329132620e-02);
    add_cartesian(sim, 5.151389e-5, 2.98721200408889800e+01, 5.18939463241125942e-01, -6.99067175203629043e-01, -2.98764379233138963e-02, 1.15313346010349438e+00, -2.30704057616862752e-02);
}

static struct reb_simulation* make_full_whfast(double dt_days) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->dt = dt_days / 365.25;
    sim->exact_finish_time = 1;
    add_full_physical_bodies(sim);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_main_solar_system(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
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
    sim->exact_finish_time = 1;
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
    sim->exact_finish_time = 1;
    reb_simulation_add_fmt(sim, "m", 9.547919e-4);
    add_cartesian(sim, 2.388e-8, 2.48120485765584205e-03, -1.34556438159649208e-03, -1.31630007347439604e-05, 1.75341752994026812e+00, 3.20163492363371871e+00, 1.39542430489110164e-01);
    add_cartesian(sim, 2.527e-8, 5.51813635549570501e-04, -4.47551321186325738e-03, -1.27610904592036212e-04, 2.85755306406471643e+00, 3.72439324002732552e-01, 7.16871253077594472e-02);
    add_cartesian(sim, 7.804e-8, 6.76397567656043287e-03, -2.29852901298470890e-03, 1.01873701583559195e-05, 7.40898140660848692e-01, 2.17447712709783314e+00, 9.37799108604124437e-02);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct bridge_state make_bridge_state(double dt_outer_days, int earth_moon_ratio, int jovian_ratio) {
    struct bridge_state state;
    memset(&state, 0, sizeof(state));
    const double dt_outer = dt_outer_days / 365.25;
    state.main_sim = make_main_solar_system();
    state.earth_moon = make_sub_earth_moon_horizons();
    state.jovian = make_sub_jovian_horizons();
    if (!state.main_sim || !state.earth_moon || !state.jovian) return state;
    state.bridge = reb_bridge_create(state.main_sim, dt_outer);
    if (!state.bridge) return state;
    reb_bridge_set_owns_main(state.bridge, 1);
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
    if (errno != 0 || end == text || *end != '\0' || value < 0 || value > INT_MAX) return -1;
    *out = (int)value;
    return 0;
}

static int run_whfast_once(double years, double dt_days, double* elapsed) {
    struct reb_simulation* sim = make_full_whfast(dt_days);
    if (!sim) return -1;
    const double start = monotonic_seconds();
    const enum REB_STATUS status = reb_simulation_integrate(sim, years);
    *elapsed = monotonic_seconds() - start;
    reb_simulation_free(sim);
    return status == REB_STATUS_SUCCESS ? 0 : -1;
}

static int run_bridge_once(double years, double dt_outer_days, int earth_moon_ratio, int jovian_ratio, double* elapsed) {
    struct bridge_state state = make_bridge_state(dt_outer_days, earth_moon_ratio, jovian_ratio);
    if (!state.bridge || !state.main_sim || !state.earth_moon || !state.jovian) {
        free_bridge_state(&state);
        return -1;
    }
    const double start = monotonic_seconds();
    const int status = reb_bridge_advance(state.bridge, years);
    *elapsed = monotonic_seconds() - start;
    free_bridge_state(&state);
    return status == 0 ? 0 : -1;
}

static void print_usage(const char* argv0) {
    fprintf(stderr, "usage: %s [--years N] [--dt-outer-days N] [--earth-moon-ratio N] [--jovian-ratio N] [--repeats N] [--warmups N]\n", argv0);
}

int main(int argc, char** argv) {
    double years = 100.0;
    double dt_outer_days = 1.0;
    int earth_moon_ratio = 5;
    int jovian_ratio = 15;
    int repeats = 3;
    int warmups = 1;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--years") == 0) {
            if (i + 1 >= argc || parse_double_value(argv[++i], &years) != 0 || years <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--dt-outer-days") == 0) {
            if (i + 1 >= argc || parse_double_value(argv[++i], &dt_outer_days) != 0 || dt_outer_days <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--earth-moon-ratio") == 0) {
            if (i + 1 >= argc || parse_int_value(argv[++i], &earth_moon_ratio) != 0 || earth_moon_ratio <= 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--jovian-ratio") == 0) {
            if (i + 1 >= argc || parse_int_value(argv[++i], &jovian_ratio) != 0 || jovian_ratio <= 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--repeats") == 0) {
            if (i + 1 >= argc || parse_int_value(argv[++i], &repeats) != 0 || repeats <= 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--warmups") == 0) {
            if (i + 1 >= argc || parse_int_value(argv[++i], &warmups) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    const double earth_moon_dt_days = dt_outer_days / (double)earth_moon_ratio;
    const double jovian_dt_days = dt_outer_days / (double)jovian_ratio;
    const double whfast_dt_days = earth_moon_dt_days < jovian_dt_days ? earth_moon_dt_days : jovian_dt_days;
    double whfast_total = 0.0;
    double bridge_total = 0.0;
    double whfast_best = 0.0;
    double bridge_best = 0.0;

    for (int i = 0; i < warmups; i++) {
        double elapsed = 0.0;
        if (run_whfast_once(years, whfast_dt_days, &elapsed) != 0) return EXIT_FAILURE;
        if (run_bridge_once(years, dt_outer_days, earth_moon_ratio, jovian_ratio, &elapsed) != 0) return EXIT_FAILURE;
    }

    for (int i = 0; i < repeats; i++) {
        double whfast_elapsed = 0.0;
        double bridge_elapsed = 0.0;
        if (run_whfast_once(years, whfast_dt_days, &whfast_elapsed) != 0) return EXIT_FAILURE;
        if (run_bridge_once(years, dt_outer_days, earth_moon_ratio, jovian_ratio, &bridge_elapsed) != 0) return EXIT_FAILURE;
        whfast_total += whfast_elapsed;
        bridge_total += bridge_elapsed;
        if (i == 0 || whfast_elapsed < whfast_best) whfast_best = whfast_elapsed;
        if (i == 0 || bridge_elapsed < bridge_best) bridge_best = bridge_elapsed;
        printf("repeat=%d whfast_seconds=%.9f bridge_seconds=%.9f speedup_whfast_over_bridge=%.9f\n", i + 1, whfast_elapsed, bridge_elapsed, bridge_elapsed > 0.0 ? whfast_elapsed / bridge_elapsed : 0.0);
    }

    const double whfast_mean = whfast_total / (double)repeats;
    const double bridge_mean = bridge_total / (double)repeats;
    printf("years=%.9f\n", years);
    printf("dt_outer_days=%.9f\n", dt_outer_days);
    printf("earth_moon_ratio=%d\n", earth_moon_ratio);
    printf("jovian_ratio=%d\n", jovian_ratio);
    printf("whfast_dt_days=%.9f\n", whfast_dt_days);
    printf("bridge_mean_seconds=%.9f\n", bridge_mean);
    printf("whfast_mean_seconds=%.9f\n", whfast_mean);
    printf("mean_speedup_whfast_over_bridge=%.9f\n", bridge_mean > 0.0 ? whfast_mean / bridge_mean : 0.0);
    printf("bridge_best_seconds=%.9f\n", bridge_best);
    printf("whfast_best_seconds=%.9f\n", whfast_best);
    printf("best_speedup_whfast_over_bridge=%.9f\n", bridge_best > 0.0 ? whfast_best / bridge_best : 0.0);
    return EXIT_SUCCESS;
}
