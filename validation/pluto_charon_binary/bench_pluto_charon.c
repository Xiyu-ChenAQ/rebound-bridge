/* SPDX-License-Identifier: GPL-3.0-or-later */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rebound_bridge.h"
#include "../common/val_cli.h"
#include "../validation_timer.h"


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

enum pc_case_kind {
    PC_CASE_MINIMAL = 0,
    PC_CASE_OUTER = 1
};

enum pc_main_index {
    PC_MAIN_PC_COM_MINIMAL = 1,
    PC_MAIN_PC_COM_OUTER = 5
};

struct circular_body {
    const char* name;
    double mass;
    double a_au;
    double phase_rad;
};

struct pc_bridge_state {
    struct reb_simulation* main_sim;
    struct reb_simulation* binary_sub;
    struct reb_bridge* bridge;
    int host_index;
};

static const struct circular_body OUTER_BODIES[] = {
    { "jupiter", 9.547919e-4, 5.2044, 0.3 },
    { "saturn", 2.8588567e-4, 9.5826, 1.8 },
    { "uranus", 4.366244e-5, 19.2184, 3.1 },
    { "neptune", 5.151389e-5, 30.1104, 4.7 }
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

static void add_pluto_charon_physical_global(struct reb_simulation* sim) {
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

static void add_circular_heliocentric_body(struct reb_simulation* sim, const struct circular_body* body) {
    const double c = cos(body->phase_rad);
    const double s = sin(body->phase_rad);
    const double v = sqrt(REB_BRIDGE_G_AU_YR_MSUN * (PC_M_SUN + body->mass) / body->a_au);
    add_cartesian(
        sim,
        body->mass,
        body->a_au * c,
        body->a_au * s,
        0.0,
        -v * s,
        v * c,
        0.0
    );
}

static void add_outer_planets(struct reb_simulation* sim) {
    const int n = (int)(sizeof(OUTER_BODIES) / sizeof(OUTER_BODIES[0]));
    for (int i = 0; i < n; i++) {
        add_circular_heliocentric_body(sim, &OUTER_BODIES[i]);
    }
}

static struct reb_simulation* make_full_whfast(enum pc_case_kind case_kind, double dt_days) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->dt = dt_days / 365.25;
    sim->exact_finish_time = 1;
    reb_simulation_add_fmt(sim, "m", PC_M_SUN);
    if (case_kind == PC_CASE_OUTER) add_outer_planets(sim);
    add_pluto_charon_physical_global(sim);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_bridge_main(enum pc_case_kind case_kind, int* host_index) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    reb_simulation_add_fmt(sim, "m", PC_M_SUN);
    if (case_kind == PC_CASE_OUTER) {
        add_outer_planets(sim);
        *host_index = PC_MAIN_PC_COM_OUTER;
    } else {
        *host_index = PC_MAIN_PC_COM_MINIMAL;
    }
    add_cartesian(sim, PC_M_BINARY, PC_R_COM_AU, 0.0, 0.0, 0.0, PC_V_COM_AU_YR, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_binary_subsystem(void) {
    struct reb_simulation* sim = reb_simulation_create();
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

static struct pc_bridge_state make_bridge_state(enum pc_case_kind case_kind, double dt_outer_days, int binary_ratio) {
    struct pc_bridge_state state;
    memset(&state, 0, sizeof(state));
    state.host_index = -1;
    state.main_sim = make_bridge_main(case_kind, &state.host_index);
    state.binary_sub = make_binary_subsystem();
    if (!state.main_sim || !state.binary_sub) return state;

    const double dt_outer = dt_outer_days / 365.25;
    state.bridge = reb_bridge_create(state.main_sim, dt_outer);
    if (!state.bridge) return state;
    reb_bridge_set_owns_main(state.bridge, 1);
    if (reb_bridge_add_subsystem(state.bridge, state.host_index, state.binary_sub, dt_outer / (double)binary_ratio, 1) != 0) {
        return state;
    }
    return state;
}

static void free_bridge_state(struct pc_bridge_state* state) {
    if (!state) return;
    if (state->bridge) {
        reb_bridge_free(state->bridge);
    } else {
        if (state->main_sim) reb_simulation_free(state->main_sim);
        if (state->binary_sub) reb_simulation_free(state->binary_sub);
    }
    memset(state, 0, sizeof(*state));
}


static int parse_case_kind(const char* text, enum pc_case_kind* out) {
    if (strcmp(text, "minimal") == 0) {
        *out = PC_CASE_MINIMAL;
        return 0;
    }
    if (strcmp(text, "outer") == 0) {
        *out = PC_CASE_OUTER;
        return 0;
    }
    return -1;
}


static const char* case_name(enum pc_case_kind case_kind) {
    switch (case_kind) {
        case PC_CASE_MINIMAL:
            return "minimal";
        case PC_CASE_OUTER:
            return "outer";
    }
    return "unknown";
}


static int run_whfast_once(enum pc_case_kind case_kind, double years, double dt_days, double* elapsed) {
    struct reb_simulation* sim = make_full_whfast(case_kind, dt_days);
    if (!sim) return -1;
    const double start = monotonic_seconds();
    const enum REB_STATUS status = reb_simulation_integrate(sim, years);
    *elapsed = monotonic_seconds() - start;
    reb_simulation_free(sim);
    return status == REB_STATUS_SUCCESS ? 0 : -1;
}

static int run_bridge_once(
    enum pc_case_kind case_kind,
    double years,
    double dt_outer_days,
    int binary_ratio,
    enum reb_bridge_method method,
    double* elapsed
) {
    struct pc_bridge_state state = make_bridge_state(case_kind, dt_outer_days, binary_ratio);
    if (!state.bridge || !state.main_sim || !state.binary_sub) {
        free_bridge_state(&state);
        return -1;
    }
    if (reb_bridge_set_method(state.bridge, method) != 0) {
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
        "usage: %s [--case minimal|outer] [--years N] [--dt-outer-days N] [--binary-ratio N] [--bridge-method order2|eos-lf4-2|eos-lf8-6-4|forward-fg4] [--repeats N] [--warmups N]\n",
        argv0
    );
}

int main(int argc, char** argv) {
    enum pc_case_kind case_kind = PC_CASE_MINIMAL;
    enum reb_bridge_method method = REB_BRIDGE_METHOD_ORDER2;
    double years = 100.0;
    double dt_outer_days = 4.0;
    int binary_ratio = 64;
    int repeats = 3;
    int warmups = 1;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--case") == 0) {
            if (i + 1 >= argc || parse_case_kind(argv[++i], &case_kind) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--years") == 0) {
            if (i + 1 >= argc || val_parse_double(argv[++i], &years) != 0 || years <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--dt-outer-days") == 0) {
            if (i + 1 >= argc || val_parse_double(argv[++i], &dt_outer_days) != 0 || dt_outer_days <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--binary-ratio") == 0) {
            if (i + 1 >= argc || val_parse_int(argv[++i], &binary_ratio) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--bridge-method") == 0) {
            if (i + 1 >= argc || val_parse_bridge_method(argv[++i], &method) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--repeats") == 0) {
            if (i + 1 >= argc || val_parse_int(argv[++i], &repeats) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--warmups") == 0) {
            if (i + 1 >= argc || val_parse_nonnegative_int(argv[++i], &warmups) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
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
        if (run_whfast_once(case_kind, years, whfast_dt_days, &elapsed) != 0) return EXIT_FAILURE;
        if (run_bridge_once(case_kind, years, dt_outer_days, binary_ratio, method, &elapsed) != 0) return EXIT_FAILURE;
    }

    for (int i = 0; i < repeats; i++) {
        double whfast_elapsed = 0.0;
        double bridge_elapsed = 0.0;
        if (run_whfast_once(case_kind, years, whfast_dt_days, &whfast_elapsed) != 0) return EXIT_FAILURE;
        if (run_bridge_once(case_kind, years, dt_outer_days, binary_ratio, method, &bridge_elapsed) != 0) return EXIT_FAILURE;
        whfast_total += whfast_elapsed;
        bridge_total += bridge_elapsed;
        if (i == 0 || whfast_elapsed < whfast_best) whfast_best = whfast_elapsed;
        if (i == 0 || bridge_elapsed < bridge_best) bridge_best = bridge_elapsed;
        printf("repeat=%d whfast_seconds=%.9f bridge_seconds=%.9f speedup_whfast_over_bridge=%.9f\n", i + 1, whfast_elapsed, bridge_elapsed, bridge_elapsed > 0.0 ? whfast_elapsed / bridge_elapsed : 0.0);
    }

    const double bridge_mean = bridge_total / (double)repeats;
    const double whfast_mean = whfast_total / (double)repeats;
    printf("case_name=%s\n", case_name(case_kind));
    printf("years=%.9f\n", years);
    printf("dt_outer_days=%.9f\n", dt_outer_days);
    printf("binary_ratio=%d\n", binary_ratio);
    printf("bridge_method=%s\n", val_bridge_method_name(method));
    printf("whfast_dt_days=%.9f\n", whfast_dt_days);
    printf("bridge_mean_seconds=%.9f\n", bridge_mean);
    printf("whfast_mean_seconds=%.9f\n", whfast_mean);
    printf("mean_speedup_whfast_over_bridge=%.9f\n", bridge_mean > 0.0 ? whfast_mean / bridge_mean : 0.0);
    printf("bridge_best_seconds=%.9f\n", bridge_best);
    printf("whfast_best_seconds=%.9f\n", whfast_best);
    printf("best_speedup_whfast_over_bridge=%.9f\n", bridge_best > 0.0 ? whfast_best / bridge_best : 0.0);
    return EXIT_SUCCESS;
}
