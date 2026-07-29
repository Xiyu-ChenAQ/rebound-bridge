/* Local-only MNRAS speed-scaling experiment. */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rebound_bridge.h"
#include "validation_timer.h"

#define PC_M_SUN 1.0
#define PC_M_PLUTO 6.552776758009927e-09
#define PC_M_CHARON 7.975981533540863e-10
#define PC_M_BINARY (PC_M_PLUTO + PC_M_CHARON)
#define PC_R_COM_AU 39.482
#define PC_V_COM_AU_YR 0.9999546351918088
#define PC_X_PLUTO_COM_AU (-1.421401166032323e-05)
#define PC_X_CHARON_COM_AU 0.00011677715758764924
#define PC_VY_PLUTO_COM_AU_YR (-0.005107263360501845)
#define PC_VY_CHARON_COM_AU_YR 0.041959420925182254

struct bridge_state {
    struct reb_simulation* main_sim;
    struct reb_simulation* binary_sub;
    struct reb_bridge* bridge;
};

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

static void add_pluto_charon_global(struct reb_simulation* sim) {
    add_cartesian(
        sim,
        PC_M_PLUTO,
        PC_R_COM_AU + PC_X_PLUTO_COM_AU,
        0.0,
        0.0,
        0.0,
        PC_V_COM_AU_YR + PC_VY_PLUTO_COM_AU_YR,
        0.0
    );
    add_cartesian(
        sim,
        PC_M_CHARON,
        PC_R_COM_AU + PC_X_CHARON_COM_AU,
        0.0,
        0.0,
        0.0,
        PC_V_COM_AU_YR + PC_VY_CHARON_COM_AU_YR,
        0.0
    );
}

static void add_synthetic_outer_bodies(struct reb_simulation* sim, int n_outer) {
    const double golden_angle = 2.39996322972865332;
    for (int i = 0; i < n_outer; i++) {
        const double a = 5.0 + 0.55 * (double)i;
        const double phase = golden_angle * (double)i + 0.17 * (double)(i % 7);
        const double c = cos(phase);
        const double s = sin(phase);
        const double mass = 1.0e-10 * (1.0 + 0.25 * (double)(i % 5));
        const double v = sqrt(REB_BRIDGE_G_AU_YR_MSUN * (PC_M_SUN + mass) / a);
        add_cartesian(sim, mass, a * c, a * s, 0.0, -v * s, v * c, 0.0);
    }
}

static struct reb_simulation* make_full_whfast(int n_outer, double dt_days) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->dt = dt_days / 365.25;
    sim->exact_finish_time = 1;
    reb_simulation_add_fmt(sim, "m", PC_M_SUN);
    add_synthetic_outer_bodies(sim, n_outer);
    add_pluto_charon_global(sim);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_bridge_main(int n_outer, int* host_index) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    reb_simulation_add_fmt(sim, "m", PC_M_SUN);
    add_synthetic_outer_bodies(sim, n_outer);
    *host_index = n_outer + 1;
    add_cartesian(sim, PC_M_BINARY, PC_R_COM_AU, 0.0, 0.0, 0.0, PC_V_COM_AU_YR, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_binary_subsystem(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    add_cartesian(sim, PC_M_PLUTO, PC_X_PLUTO_COM_AU, 0.0, 0.0, 0.0, PC_VY_PLUTO_COM_AU_YR, 0.0);
    add_cartesian(sim, PC_M_CHARON, PC_X_CHARON_COM_AU, 0.0, 0.0, 0.0, PC_VY_CHARON_COM_AU_YR, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

static void free_bridge_state(struct bridge_state* state) {
    if (!state) return;
    if (state->bridge) {
        reb_bridge_free(state->bridge);
    } else {
        if (state->main_sim) reb_simulation_free(state->main_sim);
        if (state->binary_sub) reb_simulation_free(state->binary_sub);
    }
    memset(state, 0, sizeof(*state));
}

static struct bridge_state make_bridge_state(
    int n_outer,
    double dt_outer_days,
    int binary_ratio,
    enum reb_bridge_method method
) {
    struct bridge_state state;
    int host_index = -1;
    memset(&state, 0, sizeof(state));
    state.main_sim = make_bridge_main(n_outer, &host_index);
    state.binary_sub = make_binary_subsystem();
    if (!state.main_sim || !state.binary_sub) return state;

    const double dt_outer = dt_outer_days / 365.25;
    state.bridge = reb_bridge_create(state.main_sim, dt_outer);
    if (!state.bridge) return state;
    reb_bridge_set_owns_main(state.bridge, 1);
    if (reb_bridge_set_method(state.bridge, method) != 0) return state;
    if (reb_bridge_add_subsystem(state.bridge, host_index, state.binary_sub, dt_outer / (double)binary_ratio, 1) != 0) {
        return state;
    }
    return state;
}

static int parse_int_value(const char* text, int* out) {
    char* end = NULL;
    errno = 0;
    const long value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value < 0 || value > 2147483647L) return -1;
    *out = (int)value;
    return 0;
}

static int parse_positive_int_value(const char* text, int* out) {
    if (parse_int_value(text, out) != 0 || *out <= 0) return -1;
    return 0;
}

static int parse_double_value(const char* text, double* out) {
    char* end = NULL;
    errno = 0;
    const double value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0') return -1;
    *out = value;
    return 0;
}

static int parse_bridge_method(const char* text, enum reb_bridge_method* out) {
    if (strcmp(text, "order2") == 0) {
        *out = REB_BRIDGE_METHOD_ORDER2;
        return 0;
    }
    if (strcmp(text, "eos-lf4-2") == 0) {
        *out = REB_BRIDGE_METHOD_EOS_LF4_2;
        return 0;
    }
    if (strcmp(text, "eos-lf8-6-4") == 0) {
        *out = REB_BRIDGE_METHOD_EOS_LF8_6_4;
        return 0;
    }
    return -1;
}

static const char* bridge_method_name(enum reb_bridge_method method) {
    switch (method) {
        case REB_BRIDGE_METHOD_ORDER2:
            return "order2";
        case REB_BRIDGE_METHOD_EOS_LF4_2:
            return "eos-lf4-2";
        case REB_BRIDGE_METHOD_EOS_LF8_6_4:
            return "eos-lf8-6-4";
        default:
            return "unknown";
    }
}

static int run_whfast_once(int n_outer, double years, double dt_days, double* elapsed) {
    struct reb_simulation* sim = make_full_whfast(n_outer, dt_days);
    if (!sim) return -1;
    const double start = monotonic_seconds();
    const enum REB_STATUS status = reb_simulation_integrate(sim, years);
    *elapsed = monotonic_seconds() - start;
    reb_simulation_free(sim);
    return status == REB_STATUS_SUCCESS ? 0 : -1;
}

static int run_bridge_once(
    int n_outer,
    double years,
    double dt_outer_days,
    int binary_ratio,
    enum reb_bridge_method method,
    double* elapsed
) {
    struct bridge_state state = make_bridge_state(n_outer, dt_outer_days, binary_ratio, method);
    if (!state.bridge || !state.main_sim || !state.binary_sub) {
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
    fprintf(
        stderr,
        "usage: %s [--n-outer N] [--years N] [--dt-outer-days N] [--binary-ratio N] [--bridge-method order2|eos-lf4-2|eos-lf8-6-4] [--repeats N] [--warmups N] [--csv PATH]\n",
        argv0
    );
}

int main(int argc, char** argv) {
    int n_outer = 32;
    double years = 100.0;
    double dt_outer_days = 4.0;
    int binary_ratio = 64;
    enum reb_bridge_method method = REB_BRIDGE_METHOD_EOS_LF4_2;
    int repeats = 3;
    int warmups = 1;
    const char* csv_path = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--n-outer") == 0) {
            if (i + 1 >= argc || parse_int_value(argv[++i], &n_outer) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--years") == 0) {
            if (i + 1 >= argc || parse_double_value(argv[++i], &years) != 0 || years <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--dt-outer-days") == 0) {
            if (i + 1 >= argc || parse_double_value(argv[++i], &dt_outer_days) != 0 || dt_outer_days <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--binary-ratio") == 0) {
            if (i + 1 >= argc || parse_positive_int_value(argv[++i], &binary_ratio) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--bridge-method") == 0) {
            if (i + 1 >= argc || parse_bridge_method(argv[++i], &method) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--repeats") == 0) {
            if (i + 1 >= argc || parse_positive_int_value(argv[++i], &repeats) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--warmups") == 0) {
            if (i + 1 >= argc || parse_int_value(argv[++i], &warmups) != 0) {
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

    const double whfast_dt_days = dt_outer_days / (double)binary_ratio;
    double whfast_total = 0.0;
    double bridge_total = 0.0;
    double whfast_best = 0.0;
    double bridge_best = 0.0;

    for (int i = 0; i < warmups; i++) {
        double elapsed = 0.0;
        if (run_whfast_once(n_outer, years, whfast_dt_days, &elapsed) != 0) return EXIT_FAILURE;
        if (run_bridge_once(n_outer, years, dt_outer_days, binary_ratio, method, &elapsed) != 0) return EXIT_FAILURE;
    }

    for (int i = 0; i < repeats; i++) {
        double whfast_elapsed = 0.0;
        double bridge_elapsed = 0.0;
        if (run_whfast_once(n_outer, years, whfast_dt_days, &whfast_elapsed) != 0) return EXIT_FAILURE;
        if (run_bridge_once(n_outer, years, dt_outer_days, binary_ratio, method, &bridge_elapsed) != 0) return EXIT_FAILURE;
        whfast_total += whfast_elapsed;
        bridge_total += bridge_elapsed;
        if (i == 0 || whfast_elapsed < whfast_best) whfast_best = whfast_elapsed;
        if (i == 0 || bridge_elapsed < bridge_best) bridge_best = bridge_elapsed;
        printf(
            "repeat=%d n_outer=%d whfast_seconds=%.9f bridge_seconds=%.9f speedup_whfast_over_bridge=%.9f\n",
            i + 1,
            n_outer,
            whfast_elapsed,
            bridge_elapsed,
            bridge_elapsed > 0.0 ? whfast_elapsed / bridge_elapsed : 0.0
        );
    }

    const double whfast_mean = whfast_total / (double)repeats;
    const double bridge_mean = bridge_total / (double)repeats;
    const double mean_speedup = bridge_mean > 0.0 ? whfast_mean / bridge_mean : 0.0;
    const double best_speedup = bridge_best > 0.0 ? whfast_best / bridge_best : 0.0;

    printf("n_outer=%d\n", n_outer);
    printf("years=%.9f\n", years);
    printf("dt_outer_days=%.9f\n", dt_outer_days);
    printf("binary_ratio=%d\n", binary_ratio);
    printf("bridge_method=%s\n", bridge_method_name(method));
    printf("whfast_dt_days=%.9f\n", whfast_dt_days);
    printf("whfast_mean_seconds=%.9f\n", whfast_mean);
    printf("bridge_mean_seconds=%.9f\n", bridge_mean);
    printf("mean_speedup_whfast_over_bridge=%.9f\n", mean_speedup);
    printf("whfast_best_seconds=%.9f\n", whfast_best);
    printf("bridge_best_seconds=%.9f\n", bridge_best);
    printf("best_speedup_whfast_over_bridge=%.9f\n", best_speedup);

    if (csv_path) {
        FILE* fp = fopen(csv_path, "a");
        if (!fp) {
            fprintf(stderr, "failed to open csv path: %s\n", csv_path);
            return EXIT_FAILURE;
        }
        fprintf(
            fp,
            "%d,%s,%.17g,%.17g,%d,%.17g,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            n_outer,
            bridge_method_name(method),
            years,
            dt_outer_days,
            binary_ratio,
            whfast_dt_days,
            repeats,
            warmups,
            whfast_mean,
            bridge_mean,
            mean_speedup,
            whfast_best,
            bridge_best,
            best_speedup
        );
        fclose(fp);
    }

    return EXIT_SUCCESS;
}
