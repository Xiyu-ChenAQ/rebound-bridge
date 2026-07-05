/* SPDX-License-Identifier: GPL-3.0-or-later */

#include <errno.h>
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

enum pc_bridge_advance_mode {
    PC_BRIDGE_ADVANCE_KDK = 0,
    PC_BRIDGE_ADVANCE_DKD = 1,
    PC_BRIDGE_ADVANCE_EOS_MULTIKICK = 2
};

enum pc_ias15_sync_mode {
    PC_IAS15_SYNC_OUTER = 0,
    PC_IAS15_SYNC_OUTPUT = 1
};

enum pc_main_index {
    PC_MAIN_SUN = 0,
    PC_MAIN_PC_COM_MINIMAL = 1,
    PC_MAIN_PC_COM_OUTER = 5
};

enum pc_body_index {
    PC_BODY_SUN = 0,
    PC_BODY_JUPITER = 1,
    PC_BODY_SATURN = 2,
    PC_BODY_URANUS = 3,
    PC_BODY_NEPTUNE = 4,
    PC_BODY_PLUTO = 5,
    PC_BODY_CHARON = 6,
    PC_BODY_COUNT_MINIMAL = 3,
    PC_BODY_COUNT_OUTER = 7
};

struct pc_body_state {
    double m;
    double x, y, z;
    double vx, vy, vz;
};

struct pc_pair_state {
    struct pc_body_state pluto;
    struct pc_body_state charon;
    struct pc_body_state com;
};

struct pc_bridge_state {
    struct reb_simulation* main_sim;
    struct reb_simulation* binary_sub;
    struct reb_bridge* bridge;
    int host_index;
};

struct pc_metrics {
    double total_energy;
    struct reb_vec3d total_l;
    struct reb_vec3d total_p;
    double momentum_scale;
};

struct circular_body {
    const char* name;
    double mass;
    double a_au;
    double phase_rad;
};

static const struct circular_body OUTER_BODIES[] = {
    { "jupiter", 9.547919e-4, 5.2044, 0.3 },
    { "saturn", 2.8588567e-4, 9.5826, 1.8 },
    { "uranus", 4.366244e-5, 19.2184, 3.1 },
    { "neptune", 5.151389e-5, 30.1104, 4.7 }
};

static double sqr(double x) {
    return x * x;
}

static double vec_norm3(double x, double y, double z) {
    return sqrt(x * x + y * y + z * z);
}

static double distance_state(const struct pc_body_state* a, const struct pc_body_state* b) {
    return vec_norm3(a->x - b->x, a->y - b->y, a->z - b->z);
}

static double velocity_distance_state(const struct pc_body_state* a, const struct pc_body_state* b) {
    return vec_norm3(a->vx - b->vx, a->vy - b->vy, a->vz - b->vz);
}

static double pair_distance(const struct pc_pair_state* pair) {
    return distance_state(&pair->pluto, &pair->charon);
}

static double pair_relative_error(const struct pc_pair_state* bridge, const struct pc_pair_state* ias15) {
    const double bx = bridge->charon.x - bridge->pluto.x;
    const double by = bridge->charon.y - bridge->pluto.y;
    const double bz = bridge->charon.z - bridge->pluto.z;
    const double ix = ias15->charon.x - ias15->pluto.x;
    const double iy = ias15->charon.y - ias15->pluto.y;
    const double iz = ias15->charon.z - ias15->pluto.z;
    return vec_norm3(bx - ix, by - iy, bz - iz);
}

static double pair_phase_rad(const struct pc_pair_state* pair) {
    return atan2(pair->charon.y - pair->pluto.y, pair->charon.x - pair->pluto.x);
}

static double angle_diff_rad(double a, double b) {
    return atan2(sin(a - b), cos(a - b));
}

static double rad_to_deg(double radians) {
    return radians * 180.0 / M_PI;
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

static struct reb_simulation* make_ias15_reference(enum pc_case_kind case_kind) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->exact_finish_time = 1;
    reb_simulation_add_fmt(sim, "m", PC_M_SUN);
    if (case_kind == PC_CASE_OUTER) {
        add_outer_planets(sim);
    }
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

static struct reb_simulation* make_binary_subsystem(double inner_mu_scale) {
    struct reb_simulation* sim = reb_simulation_create();
    sim->G = REB_BRIDGE_G_AU_YR_MSUN * inner_mu_scale;
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

static struct pc_body_state particle_state(const struct reb_particle* p) {
    struct pc_body_state state;
    state.m = p->m;
    state.x = p->x;
    state.y = p->y;
    state.z = p->z;
    state.vx = p->vx;
    state.vy = p->vy;
    state.vz = p->vz;
    return state;
}

static struct reb_particle state_particle(const struct pc_body_state* state) {
    struct reb_particle p;
    memset(&p, 0, sizeof(p));
    p.m = state->m;
    p.x = state->x;
    p.y = state->y;
    p.z = state->z;
    p.vx = state->vx;
    p.vy = state->vy;
    p.vz = state->vz;
    return p;
}

static struct pc_body_state pair_com(const struct pc_body_state* a, const struct pc_body_state* b) {
    const double mt = a->m + b->m;
    struct pc_body_state out;
    out.m = mt;
    out.x = (a->m * a->x + b->m * b->x) / mt;
    out.y = (a->m * a->y + b->m * b->y) / mt;
    out.z = (a->m * a->z + b->m * b->z) / mt;
    out.vx = (a->m * a->vx + b->m * b->vx) / mt;
    out.vy = (a->m * a->vy + b->m * b->vy) / mt;
    out.vz = (a->m * a->vz + b->m * b->vz) / mt;
    return out;
}

static int reconstruct_bridge_pair(const struct pc_bridge_state* state, struct pc_pair_state* out) {
    if (!state || !state->main_sim || !state->binary_sub || !out) return -1;
    const struct reb_particle* host = &state->main_sim->particles[state->host_index];
    const struct reb_particle* pluto = &state->binary_sub->particles[0];
    const struct reb_particle* charon = &state->binary_sub->particles[1];
    out->pluto = particle_state(pluto);
    out->charon = particle_state(charon);
    out->pluto.x += host->x;
    out->pluto.y += host->y;
    out->pluto.z += host->z;
    out->pluto.vx += host->vx;
    out->pluto.vy += host->vy;
    out->pluto.vz += host->vz;
    out->charon.x += host->x;
    out->charon.y += host->y;
    out->charon.z += host->z;
    out->charon.vx += host->vx;
    out->charon.vy += host->vy;
    out->charon.vz += host->vz;
    out->com = particle_state(host);
    return 0;
}

static int ias15_pluto_index(enum pc_case_kind case_kind) {
    return case_kind == PC_CASE_OUTER ? 5 : 1;
}

static int ias15_charon_index(enum pc_case_kind case_kind) {
    return case_kind == PC_CASE_OUTER ? 6 : 2;
}

static void extract_ias15_pair(const struct reb_simulation* sim, enum pc_case_kind case_kind, struct pc_pair_state* out) {
    out->pluto = particle_state(&sim->particles[ias15_pluto_index(case_kind)]);
    out->charon = particle_state(&sim->particles[ias15_charon_index(case_kind)]);
    out->com = pair_com(&out->pluto, &out->charon);
}

static int physical_body_count(enum pc_case_kind case_kind) {
    return case_kind == PC_CASE_OUTER ? PC_BODY_COUNT_OUTER : PC_BODY_COUNT_MINIMAL;
}

static void fill_bridge_bodies(const struct pc_bridge_state* state, enum pc_case_kind case_kind, struct pc_body_state* bodies) {
    struct pc_pair_state pair;
    reconstruct_bridge_pair(state, &pair);
    bodies[PC_BODY_SUN] = particle_state(&state->main_sim->particles[PC_MAIN_SUN]);
    if (case_kind == PC_CASE_OUTER) {
        bodies[PC_BODY_JUPITER] = particle_state(&state->main_sim->particles[1]);
        bodies[PC_BODY_SATURN] = particle_state(&state->main_sim->particles[2]);
        bodies[PC_BODY_URANUS] = particle_state(&state->main_sim->particles[3]);
        bodies[PC_BODY_NEPTUNE] = particle_state(&state->main_sim->particles[4]);
        bodies[PC_BODY_PLUTO] = pair.pluto;
        bodies[PC_BODY_CHARON] = pair.charon;
    } else {
        bodies[1] = pair.pluto;
        bodies[2] = pair.charon;
    }
}

static void fill_ias15_bodies(const struct reb_simulation* sim, enum pc_case_kind case_kind, struct pc_body_state* bodies) {
    const int n_bodies = physical_body_count(case_kind);
    for (int i = 0; i < n_bodies; i++) {
        bodies[i] = particle_state(&sim->particles[i]);
    }
}

static struct pc_bridge_state make_bridge_state(
    enum pc_case_kind case_kind,
    double dt_outer_days,
    int binary_ratio,
    double inner_mu_scale
) {
    struct pc_bridge_state state;
    memset(&state, 0, sizeof(state));
    state.host_index = -1;
    state.main_sim = make_bridge_main(case_kind, &state.host_index);
    state.binary_sub = make_binary_subsystem(inner_mu_scale);
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

static struct pc_metrics compute_metrics(const struct pc_body_state* bodies, int n_bodies) {
    struct pc_metrics metrics;
    memset(&metrics, 0, sizeof(metrics));

    for (int i = 0; i < n_bodies; i++) {
        const struct pc_body_state* p = &bodies[i];
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

    for (int i = 0; i < n_bodies; i++) {
        for (int j = i + 1; j < n_bodies; j++) {
            const double r = distance_state(&bodies[i], &bodies[j]);
            metrics.total_energy -= REB_BRIDGE_G_AU_YR_MSUN * bodies[i].m * bodies[j].m / r;
        }
    }
    return metrics;
}

static double relative_energy_error(double energy, double energy0) {
    const double denom = fabs(energy0) > 0.0 ? fabs(energy0) : 1.0;
    return (energy - energy0) / denom;
}

static double relative_l_error(struct reb_vec3d l, struct reb_vec3d l0) {
    const double denom = vec_norm3(l0.x, l0.y, l0.z);
    if (denom <= 0.0) return 0.0;
    return vec_norm3(l.x - l0.x, l.y - l0.y, l.z - l0.z) / denom;
}

static double relative_momentum_error(struct reb_vec3d p, struct reb_vec3d p0, double momentum_scale0) {
    const double denom = momentum_scale0 > 0.0 ? momentum_scale0 : 1.0;
    return vec_norm3(p.x - p0.x, p.y - p0.y, p.z - p0.z) / denom;
}

static struct reb_orbit binary_orbit(const struct pc_pair_state* pair) {
    return reb_orbit_from_particle(
        REB_BRIDGE_G_AU_YR_MSUN,
        state_particle(&pair->charon),
        state_particle(&pair->pluto)
    );
}

static struct reb_orbit com_sun_orbit(const struct pc_body_state* sun, const struct pc_body_state* com) {
    return reb_orbit_from_particle(
        REB_BRIDGE_G_AU_YR_MSUN,
        state_particle(com),
        state_particle(sun)
    );
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


static int parse_bridge_advance_mode(const char* text, enum pc_bridge_advance_mode* out) {
    if (strcmp(text, "kdk") == 0) {
        *out = PC_BRIDGE_ADVANCE_KDK;
        return 0;
    }
    if (strcmp(text, "dkd") == 0) {
        *out = PC_BRIDGE_ADVANCE_DKD;
        return 0;
    }
    if (strcmp(text, "eos-multikick") == 0) {
        *out = PC_BRIDGE_ADVANCE_EOS_MULTIKICK;
        return 0;
    }
    return -1;
}

static int parse_ias15_sync_mode(const char* text, enum pc_ias15_sync_mode* out) {
    if (strcmp(text, "outer") == 0) {
        *out = PC_IAS15_SYNC_OUTER;
        return 0;
    }
    if (strcmp(text, "output") == 0) {
        *out = PC_IAS15_SYNC_OUTPUT;
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

static const char* bridge_advance_mode_name(enum pc_bridge_advance_mode mode) {
    switch (mode) {
        case PC_BRIDGE_ADVANCE_KDK:
            return "kdk";
        case PC_BRIDGE_ADVANCE_DKD:
            return "dkd";
        case PC_BRIDGE_ADVANCE_EOS_MULTIKICK:
            return "eos-multikick";
    }
    return "unknown";
}

static const char* ias15_sync_mode_name(enum pc_ias15_sync_mode mode) {
    switch (mode) {
        case PC_IAS15_SYNC_OUTER:
            return "outer";
        case PC_IAS15_SYNC_OUTPUT:
            return "output";
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

static void print_usage(const char* argv0) {
    fprintf(
        stderr,
        "usage: %s [--case minimal|outer] [--years N] [--samples N] [--dt-outer-days N] [--binary-ratio N] [--bridge-method order2|eos-lf4-2|eos-lf8-6-4|forward-fg4] [--bridge-cross-kick plain|force-gradient] [--force-gradient-scale N] [--bridge-corrector none|commutator-ab] [--bridge-corrector-scale N] [--bridge-advance-mode kdk|dkd|eos-multikick] [--ias15-sync-mode outer|output] [--inner-mu-scale N] [--csv PATH] [--initial-state-csv PATH] [--force-gradient-check] [--force-gradient-check-epsilon N]\n",
        argv0
    );
}

static void write_initial_row(
    FILE* fp,
    const char* body,
    const struct pc_body_state* bridge,
    const struct pc_body_state* ias15
) {
    fprintf(
        fp,
        "%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
        body,
        bridge->m,
        ias15->m,
        bridge->m - ias15->m,
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
        distance_state(bridge, ias15),
        velocity_distance_state(bridge, ias15)
    );
}

static int write_initial_state_csv(const char* path, enum pc_case_kind case_kind) {
    if (case_kind != PC_CASE_MINIMAL) {
        fprintf(stderr, "outer case is not implemented yet\n");
        return -1;
    }

    struct pc_bridge_state bridge_state = make_bridge_state(case_kind, 4.0, 64, 1.0);
    struct reb_simulation* ias15 = make_ias15_reference(case_kind);
    if (!bridge_state.bridge || !bridge_state.main_sim || !bridge_state.binary_sub || !ias15) {
        free_bridge_state(&bridge_state);
        if (ias15) reb_simulation_free(ias15);
        return -1;
    }

    struct pc_pair_state bridge_pair;
    struct pc_pair_state ias15_pair;
    if (reconstruct_bridge_pair(&bridge_state, &bridge_pair) != 0) {
        free_bridge_state(&bridge_state);
        reb_simulation_free(ias15);
        return -1;
    }
    extract_ias15_pair(ias15, case_kind, &ias15_pair);

    FILE* fp = fopen(path, "w");
    if (!fp) {
        free_bridge_state(&bridge_state);
        reb_simulation_free(ias15);
        return -1;
    }

    fprintf(fp, "body,bridge_m,ias15_m,mass_diff,bridge_x,bridge_y,bridge_z,bridge_vx,bridge_vy,bridge_vz,ias15_x,ias15_y,ias15_z,ias15_vx,ias15_vy,ias15_vz,dr_au,dv_au_per_yr\n");
    write_initial_row(fp, "pluto", &bridge_pair.pluto, &ias15_pair.pluto);
    write_initial_row(fp, "charon", &bridge_pair.charon, &ias15_pair.charon);
    write_initial_row(fp, "pc_com", &bridge_pair.com, &ias15_pair.com);
    fclose(fp);

    const double pluto_dr = distance_state(&bridge_pair.pluto, &ias15_pair.pluto);
    const double charon_dr = distance_state(&bridge_pair.charon, &ias15_pair.charon);
    const double com_dr = distance_state(&bridge_pair.com, &ias15_pair.com);
    const double pluto_dv = velocity_distance_state(&bridge_pair.pluto, &ias15_pair.pluto);
    const double charon_dv = velocity_distance_state(&bridge_pair.charon, &ias15_pair.charon);
    const double com_dv = velocity_distance_state(&bridge_pair.com, &ias15_pair.com);
    const double pluto_dm = fabs(bridge_pair.pluto.m - ias15_pair.pluto.m);
    const double charon_dm = fabs(bridge_pair.charon.m - ias15_pair.charon.m);
    const double com_dm = fabs(bridge_pair.com.m - ias15_pair.com.m);
    double max_dr = pluto_dr > charon_dr ? pluto_dr : charon_dr;
    double max_dv = pluto_dv > charon_dv ? pluto_dv : charon_dv;
    double max_dm = pluto_dm > charon_dm ? pluto_dm : charon_dm;
    if (com_dr > max_dr) max_dr = com_dr;
    if (com_dv > max_dv) max_dv = com_dv;
    if (com_dm > max_dm) max_dm = com_dm;

    printf("max_abs_mass_diff=%.17g\n", max_dm);
    printf("max_dr_au=%.17g\n", max_dr);
    printf("max_dv_au_per_yr=%.17g\n", max_dv);
    printf("max_com_dr_au=%.17g\n", com_dr);
    printf("max_com_dv_au_per_yr=%.17g\n", com_dv);

    free_bridge_state(&bridge_state);
    reb_simulation_free(ias15);
    return 0;
}

static void write_runtime_header(FILE* fp) {
    fprintf(
        fp,
        "case_name,bridge_method,bridge_cross_kick,force_gradient_scale,bridge_corrector,bridge_corrector_scale,bridge_advance_mode,ias15_sync_mode,inner_mu_scale,dt_outer_days,binary_ratio,time_yr,"
        "bridge_step_seconds,ias15_step_seconds,step_speedup,"
        "bridge_cumulative_seconds,ias15_cumulative_seconds,cumulative_speedup,"
        "bridge_energy_rel,ias15_energy_rel,energy_rel_diff,"
        "bridge_l_rel,ias15_l_rel,l_rel_diff,"
        "bridge_p_rel,ias15_p_rel,p_rel_diff,"
        "pc_com_abs_error_au,pc_com_velocity_error_au_per_yr,"
        "pluto_abs_error_au,charon_abs_error_au,pc_relative_error_au,"
        "bridge_pc_distance_au,ias15_pc_distance_au,pc_distance_diff_au,"
        "bridge_pc_a_au,ias15_pc_a_au,pc_a_diff_au,"
        "bridge_pc_e,ias15_pc_e,pc_e_diff,"
        "bridge_pc_phase_deg,ias15_pc_phase_deg,pc_phase_diff_deg,"
        "bridge_pc_true_anomaly_deg,ias15_pc_true_anomaly_deg,pc_true_anomaly_diff_deg,"
        "bridge_pc_mean_anomaly_deg,ias15_pc_mean_anomaly_deg,pc_mean_anomaly_diff_deg,"
        "pc_phase_arc_error_au,"
        "bridge_com_sun_a_au,ias15_com_sun_a_au,com_sun_a_diff_au,"
        "bridge_com_sun_e,ias15_com_sun_e,com_sun_e_diff\n"
    );
}

static void write_runtime_row(
    FILE* fp,
    enum pc_case_kind case_kind,
    enum reb_bridge_method method,
    enum reb_bridge_cross_kick_kernel cross_kick_kernel,
    double force_gradient_scale,
    enum reb_bridge_corrector bridge_corrector,
    double bridge_corrector_scale,
    enum pc_bridge_advance_mode advance_mode,
    enum pc_ias15_sync_mode ias15_sync_mode,
    double inner_mu_scale,
    double dt_outer_days,
    int binary_ratio,
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
    double pc_com_abs_error_au,
    double pc_com_velocity_error_au_per_yr,
    double pluto_abs_error_au,
    double charon_abs_error_au,
    double pc_relative_error_au,
    const struct pc_pair_state* bridge_pair,
    const struct pc_pair_state* ias15_pair,
    double bridge_pc_distance_au,
    double ias15_pc_distance_au,
    const struct reb_orbit* bridge_pc_orbit,
    const struct reb_orbit* ias15_pc_orbit,
    const struct reb_orbit* bridge_com_orbit,
    const struct reb_orbit* ias15_com_orbit
) {
#define WRITE_CSV_VALUE(value) fprintf(fp, ",%.17g", (double)(value))
    const double bridge_phase = pair_phase_rad(bridge_pair);
    const double ias15_phase = pair_phase_rad(ias15_pair);
    const double phase_diff = angle_diff_rad(bridge_phase, ias15_phase);
    const double true_anomaly_diff = angle_diff_rad(bridge_pc_orbit->f, ias15_pc_orbit->f);
    const double mean_anomaly_diff = angle_diff_rad(bridge_pc_orbit->M, ias15_pc_orbit->M);
    const double pc_phase_arc_error_au = fabs(phase_diff) * 0.5 * (bridge_pc_distance_au + ias15_pc_distance_au);

    fprintf(
        fp,
        "%s,%s,%s",
        case_name(case_kind),
        val_bridge_method_name(method),
        cross_kick_kernel_name(cross_kick_kernel)
    );
    WRITE_CSV_VALUE(force_gradient_scale);
    fprintf(fp, ",%s", bridge_corrector_name(bridge_corrector));
    WRITE_CSV_VALUE(bridge_corrector_scale);
    fprintf(fp, ",%s,%s", bridge_advance_mode_name(advance_mode), ias15_sync_mode_name(ias15_sync_mode));
    WRITE_CSV_VALUE(inner_mu_scale);
    WRITE_CSV_VALUE(dt_outer_days);
    fprintf(fp, ",%d", binary_ratio);
    WRITE_CSV_VALUE(time_yr);
    WRITE_CSV_VALUE(bridge_step_seconds);
    WRITE_CSV_VALUE(ias15_step_seconds);
    WRITE_CSV_VALUE(bridge_step_seconds > 0.0 ? ias15_step_seconds / bridge_step_seconds : 0.0);
    WRITE_CSV_VALUE(bridge_cumulative_seconds);
    WRITE_CSV_VALUE(ias15_cumulative_seconds);
    WRITE_CSV_VALUE(bridge_cumulative_seconds > 0.0 ? ias15_cumulative_seconds / bridge_cumulative_seconds : 0.0);
    WRITE_CSV_VALUE(bridge_energy_rel);
    WRITE_CSV_VALUE(ias15_energy_rel);
    WRITE_CSV_VALUE(bridge_energy_rel - ias15_energy_rel);
    WRITE_CSV_VALUE(bridge_l_rel);
    WRITE_CSV_VALUE(ias15_l_rel);
    WRITE_CSV_VALUE(bridge_l_rel - ias15_l_rel);
    WRITE_CSV_VALUE(bridge_p_rel);
    WRITE_CSV_VALUE(ias15_p_rel);
    WRITE_CSV_VALUE(bridge_p_rel - ias15_p_rel);
    WRITE_CSV_VALUE(pc_com_abs_error_au);
    WRITE_CSV_VALUE(pc_com_velocity_error_au_per_yr);
    WRITE_CSV_VALUE(pluto_abs_error_au);
    WRITE_CSV_VALUE(charon_abs_error_au);
    WRITE_CSV_VALUE(pc_relative_error_au);
    WRITE_CSV_VALUE(bridge_pc_distance_au);
    WRITE_CSV_VALUE(ias15_pc_distance_au);
    WRITE_CSV_VALUE(bridge_pc_distance_au - ias15_pc_distance_au);
    WRITE_CSV_VALUE(bridge_pc_orbit->a);
    WRITE_CSV_VALUE(ias15_pc_orbit->a);
    WRITE_CSV_VALUE(bridge_pc_orbit->a - ias15_pc_orbit->a);
    WRITE_CSV_VALUE(bridge_pc_orbit->e);
    WRITE_CSV_VALUE(ias15_pc_orbit->e);
    WRITE_CSV_VALUE(bridge_pc_orbit->e - ias15_pc_orbit->e);
    WRITE_CSV_VALUE(rad_to_deg(bridge_phase));
    WRITE_CSV_VALUE(rad_to_deg(ias15_phase));
    WRITE_CSV_VALUE(rad_to_deg(phase_diff));
    WRITE_CSV_VALUE(rad_to_deg(bridge_pc_orbit->f));
    WRITE_CSV_VALUE(rad_to_deg(ias15_pc_orbit->f));
    WRITE_CSV_VALUE(rad_to_deg(true_anomaly_diff));
    WRITE_CSV_VALUE(rad_to_deg(bridge_pc_orbit->M));
    WRITE_CSV_VALUE(rad_to_deg(ias15_pc_orbit->M));
    WRITE_CSV_VALUE(rad_to_deg(mean_anomaly_diff));
    WRITE_CSV_VALUE(pc_phase_arc_error_au);
    WRITE_CSV_VALUE(bridge_com_orbit->a);
    WRITE_CSV_VALUE(ias15_com_orbit->a);
    WRITE_CSV_VALUE(bridge_com_orbit->a - ias15_com_orbit->a);
    WRITE_CSV_VALUE(bridge_com_orbit->e);
    WRITE_CSV_VALUE(ias15_com_orbit->e);
    WRITE_CSV_VALUE(bridge_com_orbit->e - ias15_com_orbit->e);
    fprintf(fp, "\n");
#undef WRITE_CSV_VALUE
}

static int pc_bridge_drift_all_to_target(struct pc_bridge_state* state, double target) {
    enum REB_STATUS sub_status = reb_simulation_integrate(state->binary_sub, target);
    if (sub_status != REB_STATUS_SUCCESS) return -1;

    enum REB_STATUS main_status = reb_simulation_integrate(state->main_sim, target);
    if (main_status != REB_STATUS_SUCCESS) return -1;
    return 0;
}

static double pc_time_advance_tolerance(double a, double b) {
    const double scale = fmax(fabs(a), fabs(b));
    return 1.0e-14 * fmax(1.0, scale);
}

static int pc_reached_time(double t, double target) {
    return t >= target - pc_time_advance_tolerance(t, target);
}

static int pc_bridge_advance_dkd(struct pc_bridge_state* state, double duration, double dt_outer_days) {
    if (!state || !state->bridge || !state->main_sim || !state->binary_sub) return -1;
    if (duration < 0.0) return -1;
    if (duration == 0.0) return 0;

    const double dt_outer = dt_outer_days / 365.25;
    const double target = state->main_sim->t + duration;
    while (!pc_reached_time(state->main_sim->t, target)) {
        double dt = target - state->main_sim->t;
        if (dt > dt_outer) dt = dt_outer;

        if (pc_bridge_drift_all_to_target(state, state->main_sim->t + 0.5 * dt) != 0) return -1;
        if (reb_bridge_apply_cross_kick(state->bridge, dt) != 0) return -1;
        if (pc_bridge_drift_all_to_target(state, state->main_sim->t + 0.5 * dt) != 0) return -1;
    }
    return 0;
}

static int pc_bridge_eos_multikick_step(struct pc_bridge_state* state, double dt) {
    const double a = 0.211324865405187117745425609749;
    const double mid_half = 0.5 * (1.0 - 2.0 * a);

    if (pc_bridge_drift_all_to_target(state, state->main_sim->t + a * dt) != 0) return -1;
    if (reb_bridge_apply_cross_kick(state->bridge, 0.25 * dt) != 0) return -1;
    if (pc_bridge_drift_all_to_target(state, state->main_sim->t + mid_half * dt) != 0) return -1;
    if (reb_bridge_apply_cross_kick(state->bridge, 0.5 * dt) != 0) return -1;
    if (pc_bridge_drift_all_to_target(state, state->main_sim->t + mid_half * dt) != 0) return -1;
    if (reb_bridge_apply_cross_kick(state->bridge, 0.25 * dt) != 0) return -1;
    if (pc_bridge_drift_all_to_target(state, state->main_sim->t + a * dt) != 0) return -1;
    return 0;
}

static int pc_bridge_advance_eos_multikick(struct pc_bridge_state* state, double duration, double dt_outer_days) {
    if (!state || !state->bridge || !state->main_sim || !state->binary_sub) return -1;
    if (duration < 0.0) return -1;
    if (duration == 0.0) return 0;

    const double dt_outer = dt_outer_days / 365.25;
    const double target = state->main_sim->t + duration;
    while (!pc_reached_time(state->main_sim->t, target)) {
        double dt = target - state->main_sim->t;
        if (dt > dt_outer) dt = dt_outer;
        if (pc_bridge_eos_multikick_step(state, dt) != 0) return -1;
    }
    return 0;
}

static int pc_bridge_advance(
    struct pc_bridge_state* state,
    double duration,
    double dt_outer_days,
    enum reb_bridge_method method,
    enum pc_bridge_advance_mode advance_mode
) {
    switch (advance_mode) {
        case PC_BRIDGE_ADVANCE_KDK:
            return reb_bridge_advance(state->bridge, duration);
        case PC_BRIDGE_ADVANCE_DKD:
            if (method != REB_BRIDGE_METHOD_ORDER2) {
                fprintf(stderr, "DKD advance mode is only defined for bridge-method order2\n");
                return -1;
            }
            return pc_bridge_advance_dkd(state, duration, dt_outer_days);
        case PC_BRIDGE_ADVANCE_EOS_MULTIKICK:
            if (method != REB_BRIDGE_METHOD_EOS_LF4_2) {
                fprintf(stderr, "EOS multi-kick advance mode is only defined for bridge-method eos-lf4-2\n");
                return -1;
            }
            return pc_bridge_advance_eos_multikick(state, duration, dt_outer_days);
    }
    return -1;
}

static int run_comparison(
    const char* csv_path,
    enum pc_case_kind case_kind,
    double years,
    int samples,
    double dt_outer_days,
    int binary_ratio,
    double inner_mu_scale,
    enum pc_ias15_sync_mode ias15_sync_mode,
    enum reb_bridge_method method,
    enum reb_bridge_cross_kick_kernel cross_kick_kernel,
    double force_gradient_scale,
    enum reb_bridge_corrector bridge_corrector,
    double bridge_corrector_scale,
    enum pc_bridge_advance_mode advance_mode
) {
    struct pc_bridge_state bridge_state = make_bridge_state(case_kind, dt_outer_days, binary_ratio, inner_mu_scale);
    struct reb_simulation* ias15 = make_ias15_reference(case_kind);
    if (!bridge_state.bridge || !bridge_state.main_sim || !bridge_state.binary_sub || !ias15) {
        free_bridge_state(&bridge_state);
        if (ias15) reb_simulation_free(ias15);
        return -1;
    }
    if (reb_bridge_set_method(bridge_state.bridge, method) != 0) {
        free_bridge_state(&bridge_state);
        reb_simulation_free(ias15);
        return -1;
    }
    if (reb_bridge_set_cross_kick_kernel(bridge_state.bridge, cross_kick_kernel) != 0) {
        free_bridge_state(&bridge_state);
        reb_simulation_free(ias15);
        return -1;
    }
    if (reb_bridge_set_force_gradient_scale(bridge_state.bridge, force_gradient_scale) != 0) {
        free_bridge_state(&bridge_state);
        reb_simulation_free(ias15);
        return -1;
    }
    if (reb_bridge_set_corrector(bridge_state.bridge, bridge_corrector) != 0 ||
        reb_bridge_set_corrector_scale(bridge_state.bridge, bridge_corrector_scale) != 0) {
        free_bridge_state(&bridge_state);
        reb_simulation_free(ias15);
        return -1;
    }

    FILE* fp = fopen(csv_path, "w");
    if (!fp) {
        free_bridge_state(&bridge_state);
        reb_simulation_free(ias15);
        return -1;
    }
    write_runtime_header(fp);

    const int n_bodies = physical_body_count(case_kind);
    struct pc_body_state bridge_bodies0[PC_BODY_COUNT_OUTER];
    struct pc_body_state ias15_bodies0[PC_BODY_COUNT_OUTER];
    fill_bridge_bodies(&bridge_state, case_kind, bridge_bodies0);
    fill_ias15_bodies(ias15, case_kind, ias15_bodies0);
    const struct pc_metrics bridge_metrics0 = compute_metrics(bridge_bodies0, n_bodies);
    const struct pc_metrics ias15_metrics0 = compute_metrics(ias15_bodies0, n_bodies);

    if (bridge_corrector != REB_BRIDGE_CORRECTOR_NONE &&
        reb_bridge_apply_corrector(bridge_state.bridge, REB_BRIDGE_CORRECTOR_FORWARD) != 0) {
        fclose(fp);
        free_bridge_state(&bridge_state);
        reb_simulation_free(ias15);
        return -1;
    }

    double bridge_total_seconds = 0.0;
    double ias15_total_seconds = 0.0;
    double max_pc_com_abs_error_au = 0.0;
    double max_pc_com_velocity_error_au_per_yr = 0.0;
    double max_pluto_abs_error_au = 0.0;
    double max_charon_abs_error_au = 0.0;
    double max_pc_relative_error_au = 0.0;
    double max_abs_pc_distance_diff_au = 0.0;
    double max_abs_pc_a_diff_au = 0.0;
    double max_abs_pc_e_diff = 0.0;
    double max_abs_pc_phase_arc_error_au = 0.0;
    double max_abs_bridge_energy_rel = 0.0;
    double max_abs_bridge_l_rel = 0.0;
    double max_abs_bridge_p_rel = 0.0;

    const double dt_outer = dt_outer_days / 365.25;
    const int total_events = (int)ceil(years / dt_outer - 1.0e-14);
    int next_sample = 1;
    double bridge_seconds_since_output = 0.0;
    double ias15_seconds_since_output = 0.0;

    for (int event = 1; event <= total_events; event++) {
        double step_duration = years - bridge_state.main_sim->t;
        if (step_duration > dt_outer) step_duration = dt_outer;
        if (step_duration <= 0.0) break;

        const double bridge_start = monotonic_seconds();
        if (pc_bridge_advance(&bridge_state, step_duration, dt_outer_days, method, advance_mode) != 0) {
            fclose(fp);
            free_bridge_state(&bridge_state);
            reb_simulation_free(ias15);
            return -1;
        }
        const double bridge_event_seconds = monotonic_seconds() - bridge_start;
        bridge_total_seconds += bridge_event_seconds;
        bridge_seconds_since_output += bridge_event_seconds;

        if (ias15_sync_mode == PC_IAS15_SYNC_OUTER) {
            const double ias15_start = monotonic_seconds();
            const enum REB_STATUS ias15_status = reb_simulation_integrate(ias15, bridge_state.main_sim->t);
            const double ias15_event_seconds = monotonic_seconds() - ias15_start;
            if (ias15_status != REB_STATUS_SUCCESS) {
                fclose(fp);
                free_bridge_state(&bridge_state);
                reb_simulation_free(ias15);
                return -1;
            }
            ias15_total_seconds += ias15_event_seconds;
            ias15_seconds_since_output += ias15_event_seconds;
        }

        if ((long long)event * (long long)samples < (long long)next_sample * (long long)total_events) {
            continue;
        }

        const double output_time = bridge_state.main_sim->t;
        const double bridge_step_seconds = bridge_seconds_since_output;
        bridge_seconds_since_output = 0.0;

        if (ias15_sync_mode == PC_IAS15_SYNC_OUTPUT) {
            const double ias15_start = monotonic_seconds();
            const enum REB_STATUS ias15_status = reb_simulation_integrate(ias15, output_time);
            const double ias15_output_seconds = monotonic_seconds() - ias15_start;
            if (ias15_status != REB_STATUS_SUCCESS) {
                fclose(fp);
                free_bridge_state(&bridge_state);
                reb_simulation_free(ias15);
                return -1;
            }
            ias15_total_seconds += ias15_output_seconds;
            ias15_seconds_since_output += ias15_output_seconds;
        }

        const double ias15_step_seconds = ias15_seconds_since_output;
        ias15_seconds_since_output = 0.0;

        if (bridge_corrector != REB_BRIDGE_CORRECTOR_NONE &&
            reb_bridge_apply_corrector(bridge_state.bridge, REB_BRIDGE_CORRECTOR_INVERSE) != 0) {
            fclose(fp);
            free_bridge_state(&bridge_state);
            reb_simulation_free(ias15);
            return -1;
        }

        struct pc_pair_state bridge_pair;
        struct pc_pair_state ias15_pair;
        struct pc_body_state bridge_bodies[PC_BODY_COUNT_OUTER];
        struct pc_body_state ias15_bodies[PC_BODY_COUNT_OUTER];
        reconstruct_bridge_pair(&bridge_state, &bridge_pair);
        extract_ias15_pair(ias15, case_kind, &ias15_pair);
        fill_bridge_bodies(&bridge_state, case_kind, bridge_bodies);
        fill_ias15_bodies(ias15, case_kind, ias15_bodies);

        const struct pc_metrics bridge_metrics = compute_metrics(bridge_bodies, n_bodies);
        const struct pc_metrics ias15_metrics = compute_metrics(ias15_bodies, n_bodies);
        const double bridge_energy_rel = relative_energy_error(bridge_metrics.total_energy, bridge_metrics0.total_energy);
        const double ias15_energy_rel = relative_energy_error(ias15_metrics.total_energy, ias15_metrics0.total_energy);
        const double bridge_l_rel = relative_l_error(bridge_metrics.total_l, bridge_metrics0.total_l);
        const double ias15_l_rel = relative_l_error(ias15_metrics.total_l, ias15_metrics0.total_l);
        const double bridge_p_rel = relative_momentum_error(bridge_metrics.total_p, bridge_metrics0.total_p, bridge_metrics0.momentum_scale);
        const double ias15_p_rel = relative_momentum_error(ias15_metrics.total_p, ias15_metrics0.total_p, ias15_metrics0.momentum_scale);

        const double pc_com_abs_error_au = distance_state(&bridge_pair.com, &ias15_pair.com);
        const double pc_com_velocity_error_au_per_yr = velocity_distance_state(&bridge_pair.com, &ias15_pair.com);
        const double pluto_abs_error_au = distance_state(&bridge_pair.pluto, &ias15_pair.pluto);
        const double charon_abs_error_au = distance_state(&bridge_pair.charon, &ias15_pair.charon);
        const double pc_relative_error_au = pair_relative_error(&bridge_pair, &ias15_pair);
        const double bridge_pc_distance_au = pair_distance(&bridge_pair);
        const double ias15_pc_distance_au = pair_distance(&ias15_pair);
        const struct reb_orbit bridge_pc_orbit = binary_orbit(&bridge_pair);
        const struct reb_orbit ias15_pc_orbit = binary_orbit(&ias15_pair);
        const struct reb_orbit bridge_com_orbit = com_sun_orbit(&bridge_bodies[PC_BODY_SUN], &bridge_pair.com);
        const struct reb_orbit ias15_com_orbit = com_sun_orbit(&ias15_bodies[PC_BODY_SUN], &ias15_pair.com);
        const double pc_phase_arc_error_au =
            fabs(angle_diff_rad(pair_phase_rad(&bridge_pair), pair_phase_rad(&ias15_pair))) *
            0.5 * (bridge_pc_distance_au + ias15_pc_distance_au);

        if (pc_com_abs_error_au > max_pc_com_abs_error_au) max_pc_com_abs_error_au = pc_com_abs_error_au;
        if (pc_com_velocity_error_au_per_yr > max_pc_com_velocity_error_au_per_yr) max_pc_com_velocity_error_au_per_yr = pc_com_velocity_error_au_per_yr;
        if (pluto_abs_error_au > max_pluto_abs_error_au) max_pluto_abs_error_au = pluto_abs_error_au;
        if (charon_abs_error_au > max_charon_abs_error_au) max_charon_abs_error_au = charon_abs_error_au;
        if (pc_relative_error_au > max_pc_relative_error_au) max_pc_relative_error_au = pc_relative_error_au;
        if (fabs(bridge_pc_distance_au - ias15_pc_distance_au) > max_abs_pc_distance_diff_au) max_abs_pc_distance_diff_au = fabs(bridge_pc_distance_au - ias15_pc_distance_au);
        if (fabs(bridge_pc_orbit.a - ias15_pc_orbit.a) > max_abs_pc_a_diff_au) max_abs_pc_a_diff_au = fabs(bridge_pc_orbit.a - ias15_pc_orbit.a);
        if (fabs(bridge_pc_orbit.e - ias15_pc_orbit.e) > max_abs_pc_e_diff) max_abs_pc_e_diff = fabs(bridge_pc_orbit.e - ias15_pc_orbit.e);
        if (pc_phase_arc_error_au > max_abs_pc_phase_arc_error_au) max_abs_pc_phase_arc_error_au = pc_phase_arc_error_au;
        if (fabs(bridge_energy_rel) > max_abs_bridge_energy_rel) max_abs_bridge_energy_rel = fabs(bridge_energy_rel);
        if (fabs(bridge_l_rel) > max_abs_bridge_l_rel) max_abs_bridge_l_rel = fabs(bridge_l_rel);
        if (fabs(bridge_p_rel) > max_abs_bridge_p_rel) max_abs_bridge_p_rel = fabs(bridge_p_rel);

        write_runtime_row(
            fp,
            case_kind,
            method,
            cross_kick_kernel,
            force_gradient_scale,
            bridge_corrector,
            bridge_corrector_scale,
            advance_mode,
            ias15_sync_mode,
            inner_mu_scale,
            dt_outer_days,
            binary_ratio,
            output_time,
            bridge_step_seconds,
            ias15_step_seconds,
            bridge_total_seconds,
            ias15_total_seconds,
            bridge_energy_rel,
            ias15_energy_rel,
            bridge_l_rel,
            ias15_l_rel,
            bridge_p_rel,
            ias15_p_rel,
            pc_com_abs_error_au,
            pc_com_velocity_error_au_per_yr,
            pluto_abs_error_au,
            charon_abs_error_au,
            pc_relative_error_au,
            &bridge_pair,
            &ias15_pair,
            bridge_pc_distance_au,
            ias15_pc_distance_au,
            &bridge_pc_orbit,
            &ias15_pc_orbit,
            &bridge_com_orbit,
            &ias15_com_orbit
        );

        if (bridge_corrector != REB_BRIDGE_CORRECTOR_NONE &&
            reb_bridge_apply_corrector(bridge_state.bridge, REB_BRIDGE_CORRECTOR_FORWARD) != 0) {
            fclose(fp);
            free_bridge_state(&bridge_state);
            reb_simulation_free(ias15);
            return -1;
        }
        next_sample++;
    }

    fclose(fp);

    printf("case_name=%s\n", case_name(case_kind));
    printf("bridge_method=%s\n", val_bridge_method_name(method));
    printf("bridge_cross_kick=%s\n", cross_kick_kernel_name(cross_kick_kernel));
    printf("force_gradient_scale=%.17g\n", force_gradient_scale);
    printf("bridge_corrector=%s\n", bridge_corrector_name(bridge_corrector));
    printf("bridge_corrector_scale=%.17g\n", bridge_corrector_scale);
    printf("bridge_advance_mode=%s\n", bridge_advance_mode_name(advance_mode));
    printf("ias15_sync_mode=%s\n", ias15_sync_mode_name(ias15_sync_mode));
    printf("inner_mu_scale=%.17g\n", inner_mu_scale);
    printf("years=%.17g\n", years);
    printf("samples=%d\n", samples);
    printf("dt_outer_days=%.17g\n", dt_outer_days);
    printf("binary_ratio=%d\n", binary_ratio);
    printf("binary_dt_days=%.17g\n", dt_outer_days / (double)binary_ratio);
    printf("bridge_total_seconds=%.17g\n", bridge_total_seconds);
    printf("ias15_total_seconds=%.17g\n", ias15_total_seconds);
    printf("overall_speedup=%.17g\n", bridge_total_seconds > 0.0 ? ias15_total_seconds / bridge_total_seconds : 0.0);
    printf("max_pc_com_abs_error_au=%.17g\n", max_pc_com_abs_error_au);
    printf("max_pc_com_velocity_error_au_per_yr=%.17g\n", max_pc_com_velocity_error_au_per_yr);
    printf("max_pluto_abs_error_au=%.17g\n", max_pluto_abs_error_au);
    printf("max_charon_abs_error_au=%.17g\n", max_charon_abs_error_au);
    printf("max_pc_relative_error_au=%.17g\n", max_pc_relative_error_au);
    printf("max_abs_pc_distance_diff_au=%.17g\n", max_abs_pc_distance_diff_au);
    printf("max_abs_pc_a_diff_au=%.17g\n", max_abs_pc_a_diff_au);
    printf("max_abs_pc_e_diff=%.17g\n", max_abs_pc_e_diff);
    printf("max_abs_pc_phase_arc_error_au=%.17g\n", max_abs_pc_phase_arc_error_au);
    printf("max_abs_bridge_energy_rel=%.17g\n", max_abs_bridge_energy_rel);
    printf("max_abs_bridge_l_rel=%.17g\n", max_abs_bridge_l_rel);
    printf("max_abs_bridge_p_rel=%.17g\n", max_abs_bridge_p_rel);

    free_bridge_state(&bridge_state);
    reb_simulation_free(ias15);
    return 0;
}

static int run_force_gradient_check(
    enum pc_case_kind case_kind,
    double dt_outer_days,
    int binary_ratio,
    double inner_mu_scale,
    enum reb_bridge_method method,
    double epsilon
) {
    struct pc_bridge_state bridge_state = make_bridge_state(case_kind, dt_outer_days, binary_ratio, inner_mu_scale);
    if (!bridge_state.bridge || !bridge_state.main_sim || !bridge_state.binary_sub) {
        free_bridge_state(&bridge_state);
        return -1;
    }
    if (reb_bridge_set_method(bridge_state.bridge, method) != 0) {
        free_bridge_state(&bridge_state);
        return -1;
    }

    struct reb_bridge_force_gradient_check check;
    if (reb_bridge_force_gradient_finite_difference_check(
            bridge_state.bridge,
            0,
            epsilon,
            &check
        ) != 0) {
        free_bridge_state(&bridge_state);
        return -1;
    }

    printf("force_gradient_check_epsilon=%.17g\n", epsilon);
    printf("force_gradient_max_abs_difference=%.17g\n", check.max_abs_difference);
    printf("force_gradient_max_relative_difference=%.17g\n", check.max_relative_difference);
    printf("force_gradient_max_abs_analytic=%.17g\n", check.max_abs_analytic);
    printf("force_gradient_max_abs_finite_difference=%.17g\n", check.max_abs_finite_difference);
    printf("force_gradient_internal_weighted_acc_residual=%.17g\n", check.force_gradient_internal_weighted_acc_residual);
    printf("force_gradient_main_weighted_acc_residual=%.17g\n", check.force_gradient_main_weighted_acc_residual);
    printf("cross_kick_internal_weighted_acc_residual=%.17g\n", check.internal_weighted_acc_residual);
    printf("cross_kick_main_weighted_acc_residual=%.17g\n", check.main_weighted_acc_residual);

    free_bridge_state(&bridge_state);
    return 0;
}

int main(int argc, char** argv) {
    enum pc_case_kind case_kind = PC_CASE_MINIMAL;
    enum reb_bridge_method method = REB_BRIDGE_METHOD_ORDER2;
    enum reb_bridge_cross_kick_kernel cross_kick_kernel = REB_BRIDGE_CROSS_KICK_PLAIN;
    enum reb_bridge_corrector bridge_corrector = REB_BRIDGE_CORRECTOR_NONE;
    enum pc_bridge_advance_mode advance_mode = PC_BRIDGE_ADVANCE_KDK;
    enum pc_ias15_sync_mode ias15_sync_mode = PC_IAS15_SYNC_OUTER;
    const char* initial_state_csv = NULL;
    const char* csv_path = NULL;
    int force_gradient_check = 0;
    double years = 100.0;
    int samples = 200;
    double dt_outer_days = 4.0;
    double force_gradient_scale = 1.0;
    double force_gradient_check_epsilon = 1.0e-7;
    double bridge_corrector_scale = 1.0;
    double inner_mu_scale = 1.0;
    int binary_ratio = 64;

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
            if (i + 1 >= argc || val_parse_double(argv[++i], &force_gradient_scale) != 0 || !isfinite(force_gradient_scale)) {
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
            if (i + 1 >= argc || val_parse_double(argv[++i], &bridge_corrector_scale) != 0 || !isfinite(bridge_corrector_scale)) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--bridge-advance-mode") == 0) {
            if (i + 1 >= argc || parse_bridge_advance_mode(argv[++i], &advance_mode) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--ias15-sync-mode") == 0) {
            if (i + 1 >= argc || parse_ias15_sync_mode(argv[++i], &ias15_sync_mode) != 0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--inner-mu-scale") == 0) {
            if (i + 1 >= argc || val_parse_double(argv[++i], &inner_mu_scale) != 0 || inner_mu_scale <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[i], "--csv") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            csv_path = argv[++i];
        } else if (strcmp(argv[i], "--initial-state-csv") == 0) {
            if (i + 1 >= argc) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            initial_state_csv = argv[++i];
        } else if (strcmp(argv[i], "--force-gradient-check") == 0) {
            force_gradient_check = 1;
        } else if (strcmp(argv[i], "--force-gradient-check-epsilon") == 0) {
            if (i + 1 >= argc || val_parse_double(argv[++i], &force_gradient_check_epsilon) != 0 || force_gradient_check_epsilon <= 0.0) {
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else {
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (initial_state_csv) {
        return write_initial_state_csv(initial_state_csv, case_kind) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    if (force_gradient_check) {
        return run_force_gradient_check(
            case_kind,
            dt_outer_days,
            binary_ratio,
            inner_mu_scale,
            method,
            force_gradient_check_epsilon
        ) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    if (csv_path) {
        return run_comparison(csv_path, case_kind, years, samples, dt_outer_days, binary_ratio, inner_mu_scale, ias15_sync_mode, method, cross_kick_kernel, force_gradient_scale, bridge_corrector, bridge_corrector_scale, advance_mode) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    print_usage(argv[0]);
    return EXIT_FAILURE;
}
