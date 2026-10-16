/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#include "bridge_internal.h"

#include <float.h>
#include <stdarg.h>

/* Core Bridge lifecycle and stable API. */

static REB_BRIDGE_THREAD_LOCAL enum reb_bridge_error reb_bridge_error_state = REB_BRIDGE_ERROR_NONE;
static REB_BRIDGE_THREAD_LOCAL char reb_bridge_error_message[512];

enum reb_bridge_error reb_bridge_last_error(void) {
    return reb_bridge_error_state;
}

const char* reb_bridge_last_error_message(void) {
    return reb_bridge_error_message;
}

void reb_bridge_clear_error(void) {
    reb_bridge_error_state = REB_BRIDGE_ERROR_NONE;
    reb_bridge_error_message[0] = '\0';
}

int bridge_error_code(enum reb_bridge_error code, const char* msg) {
    reb_bridge_error_state = code;
    snprintf(reb_bridge_error_message, sizeof(reb_bridge_error_message), "%s", msg);
    fprintf(stderr, "[rebound_bridge] %s\n", reb_bridge_error_message);
    return -1;
}

int bridge_error_format(enum reb_bridge_error code, const char* format, ...) {
    va_list args;
    va_start(args, format);
    vsnprintf(reb_bridge_error_message, sizeof(reb_bridge_error_message), format, args);
    va_end(args);
    reb_bridge_error_state = code;
    fprintf(stderr, "[rebound_bridge] %s\n", reb_bridge_error_message);
    return -1;
}

static int validate_simulation_state(const struct reb_simulation* sim, const char* name) {
    if (!isfinite(sim->G) || sim->G < 0.0) {
        return bridge_error_format(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                   "%s G must be finite and non-negative", name);
    }
    if (!isfinite(sim->t)) {
        return bridge_error_format(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                   "%s time must be finite", name);
    }
    for (int i = 0; i < sim->N; i++) {
        const struct reb_particle* p = &sim->particles[i];
        if (!isfinite(p->m) || p->m < 0.0) {
            return bridge_error_format(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                       "%s particle %d mass must be finite and non-negative", name, i);
        }
        if (!isfinite(p->x) || !isfinite(p->y) || !isfinite(p->z) ||
            !isfinite(p->vx) || !isfinite(p->vy) || !isfinite(p->vz)) {
            return bridge_error_format(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                       "%s particle %d position and velocity must be finite", name, i);
        }
    }
    return 0;
}

int reb_bridge_set_integrator_whfast(struct reb_simulation* sim) {
    reb_bridge_clear_error();
    if (!sim) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "sim is NULL");
    sim->integrator = REB_INTEGRATOR_WHFAST;
    sim->ri_whfast.safe_mode = 1;
    return 0;
}

int reb_bridge_set_integrator_ias15(struct reb_simulation* sim) {
    reb_bridge_clear_error();
    if (!sim) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "sim is NULL");
    sim->integrator = REB_INTEGRATOR_IAS15;
    return 0;
}

struct reb_bridge* reb_bridge_create(struct reb_simulation* main_sim, double dt_outer) {
    reb_bridge_clear_error();
    if (!main_sim) {
        bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "main_sim is NULL");
        return NULL;
    }
    if (!isfinite(dt_outer) || dt_outer <= 0.0) {
        bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "dt_outer must be finite and positive");
        return NULL;
    }
    if (validate_simulation_state(main_sim, "main_sim") != 0) return NULL;

    struct reb_bridge* bridge = (struct reb_bridge*)calloc(1, sizeof(*bridge));
    if (!bridge) {
        bridge_error_code(REB_BRIDGE_ERROR_OUT_OF_MEMORY, "cannot allocate Bridge object");
        return NULL;
    }

    bridge->main_sim = main_sim;
    bridge->whfast_jacobi_updates = 1;
    bridge->dt_outer = dt_outer;
    bridge->owns_main_sim = 0;
    bridge->method = REB_BRIDGE_METHOD_ORDER2;
    bridge->subsystem_coupling = REB_BRIDGE_SUBSYSTEM_COUPLING_MONOPOLE;
    bridge->main_sim->dt = dt_outer;
    return bridge;
}

void reb_bridge_set_owns_main(struct reb_bridge* bridge, int owns_main_sim) {
    reb_bridge_clear_error();
    if (!bridge) {
        bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "bridge is NULL");
        return;
    }
    bridge->owns_main_sim = owns_main_sim ? 1 : 0;
}

static int valid_bridge_method(enum reb_bridge_method method) {
    switch (method) {
        case REB_BRIDGE_METHOD_ORDER2:
        case REB_BRIDGE_METHOD_EOS_LF4_2:
            return 1;
        case REB_BRIDGE_METHOD_INVALID:
            return 0;
    }
    return 0;
}

int reb_bridge_set_method(struct reb_bridge* bridge, enum reb_bridge_method method) {
    reb_bridge_clear_error();
    if (!bridge) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "bridge is NULL");
    if (!valid_bridge_method(method)) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "invalid bridge method");
    bridge->method = method;
    return 0;
}

enum reb_bridge_method reb_bridge_get_method(const struct reb_bridge* bridge) {
    return bridge ? bridge->method : REB_BRIDGE_METHOD_INVALID;
}

static int valid_subsystem_coupling(enum reb_bridge_subsystem_coupling coupling) {
    switch (coupling) {
        case REB_BRIDGE_SUBSYSTEM_COUPLING_MONOPOLE:
        case REB_BRIDGE_SUBSYSTEM_COUPLING_RESOLVED:
            return 1;
        case REB_BRIDGE_SUBSYSTEM_COUPLING_INVALID:
            return 0;
    }
    return 0;
}

int reb_bridge_set_subsystem_coupling(
    struct reb_bridge* bridge,
    enum reb_bridge_subsystem_coupling coupling
) {
    reb_bridge_clear_error();
    if (!bridge) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "bridge is NULL");
    if (!valid_subsystem_coupling(coupling)) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "invalid subsystem coupling");
    bridge->subsystem_coupling = coupling;
    return 0;
}

enum reb_bridge_subsystem_coupling reb_bridge_get_subsystem_coupling(
    const struct reb_bridge* bridge
) {
    return bridge ? bridge->subsystem_coupling : REB_BRIDGE_SUBSYSTEM_COUPLING_INVALID;
}

int reb_bridge_set_whfast_jacobi_updates(struct reb_bridge* bridge, int enabled) {
    reb_bridge_clear_error();
    if (!bridge) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "bridge is NULL");
    if (enabled != 0 && enabled != 1) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "WHFast Jacobi updates must be 0 or 1");
    }
    if (bridge->whfast_jacobi_updates != enabled) {
        bridge->whfast_jacobi_updates = enabled;
        bridge->main_whfast_cache.valid = 0;
        for (int i = 0; i < bridge->n_subsystems; i++) {
            bridge->subsystems[i].whfast_cache.valid = 0;
        }
    }
    return 0;
}

int reb_bridge_get_whfast_jacobi_updates(const struct reb_bridge* bridge) {
    return bridge ? bridge->whfast_jacobi_updates : -1;
}

static int validate_subsystem_canonical_state(
    const struct reb_particle* host,
    const struct reb_simulation* sub_sim,
    int check_barycenter
) {
    double total_mass = 0.0;
    double position_scale = 0.0;
    double velocity_scale = 0.0;
    struct reb_bridge_vec3 weighted_position = vec3(0.0, 0.0, 0.0);
    struct reb_bridge_vec3 weighted_velocity = vec3(0.0, 0.0, 0.0);

    for (int i = 0; i < sub_sim->N; i++) {
        const struct reb_particle* particle = &sub_sim->particles[i];
        if (!isfinite(particle->m) || particle->m < 0.0) {
            return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "subsystem particle masses must be finite and non-negative");
        }
        const struct reb_bridge_vec3 position = particle_pos(particle);
        const struct reb_bridge_vec3 velocity = vec3(particle->vx, particle->vy, particle->vz);
        total_mass += particle->m;
        weighted_position = vec_add(weighted_position, vec_scale(position, particle->m));
        weighted_velocity = vec_add(weighted_velocity, vec_scale(velocity, particle->m));
        position_scale = fmax(position_scale, vec_norm(position));
        velocity_scale = fmax(velocity_scale, vec_norm(velocity));
    }

    if (!isfinite(total_mass) || total_mass <= 0.0) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "subsystem total mass must be positive and finite");
    }
    if (!isfinite(host->m) || host->m <= 0.0) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "subsystem host mass must be positive and finite");
    }

    const double summation_factor = 128.0 * DBL_EPSILON * fmax(1.0, (double)sub_sim->N);
    const double mass_tolerance = summation_factor * fmax(fabs(host->m), total_mass);
    if (fabs(host->m - total_mass) > mass_tolerance) {
        return bridge_error_format(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                   "subsystem host mass %.17g must equal subsystem total mass %.17g",
                                   host->m, total_mass);
    }

    const double position_tolerance = summation_factor * total_mass * position_scale;
    if (!isfinite(vec_norm(weighted_position)) ||
        (check_barycenter && vec_norm(weighted_position) > position_tolerance)) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "subsystem barycenter must be finite and at the origin; use reb_simulation_move_to_com");
    }

    const double velocity_tolerance = summation_factor * total_mass * velocity_scale;
    if (!isfinite(vec_norm(weighted_velocity)) ||
        (check_barycenter && vec_norm(weighted_velocity) > velocity_tolerance)) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "subsystem center-of-mass velocity must be finite and zero; use reb_simulation_move_to_com");
    }
    return 0;
}

void reb_bridge_free(struct reb_bridge* bridge) {
    if (!bridge) return;
    for (int i = 0; i < bridge->n_subsystems; i++) {
        free(bridge->subsystems[i].whfast_cache.state);
        if (bridge->subsystems[i].owns_sub_sim && bridge->subsystems[i].sub_sim) {
            reb_simulation_free(bridge->subsystems[i].sub_sim);
        }
    }
    free(bridge->subsystems);
    free(bridge->main_whfast_cache.state);
    free(bridge->internal_acc);
    free(bridge->source_acc);
    free(bridge->main_acc);
    free(bridge->fg_internal_acc);
    free(bridge->fg_main_acc);
    free(bridge->fd_internal_acc);
    free(bridge->fd_main_acc);
    free(bridge->saved_sub_pos);
    free(bridge->saved_main_pos);
    if (bridge->owns_main_sim && bridge->main_sim) {
        reb_simulation_free(bridge->main_sim);
    }
    free(bridge);
}

int reb_bridge_internal_add_subsystem(
    struct reb_bridge* bridge,
    int host_index,
    struct reb_simulation* sub_sim,
    double dt_inner,
    int owns_sub_sim,
    int check_barycenter
) {
    reb_bridge_clear_error();
    if (!bridge) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "bridge is NULL");
    if (!sub_sim) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "sub_sim is NULL");
    if (sub_sim == bridge->main_sim) {
        return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                 "sub_sim must be a different simulation from main_sim");
    }
    if (!valid_particle_index(bridge->main_sim, host_index)) {
        return bridge_error_format(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                   "host_index %d is outside main_sim particle range [0, %d)",
                                   host_index, bridge->main_sim->N);
    }
    if (sub_sim->N <= 0) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "sub_sim has no particles");
    for (int i = 0; i < bridge->n_subsystems; i++) {
        if (bridge->subsystems[i].sub_sim == sub_sim ||
            bridge->subsystems[i].host_index == host_index) {
            return bridge_error_code(REB_BRIDGE_ERROR_DUPLICATE_REGISTRATION,
                                     "duplicate subsystem registration: each simulation and host_index may be registered once");
        }
    }
    if (validate_simulation_state(bridge->main_sim, "main_sim") != 0 ||
        validate_simulation_state(sub_sim, "sub_sim") != 0) return -1;
    const double G_scale = fmax(fabs(sub_sim->G), fabs(bridge->main_sim->G));
    if (fabs(sub_sim->G - bridge->main_sim->G) > 64.0 * DBL_EPSILON * G_scale) {
        return bridge_error_format(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                   "sub_sim G %.17g must match main_sim G %.17g; use consistent units",
                                   sub_sim->G, bridge->main_sim->G);
    }
    if (fabs(sub_sim->t - bridge->main_sim->t) > bridge_time_sync_tolerance(sub_sim->t, bridge->main_sim->t)) {
        return bridge_error_format(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                   "sub_sim time %.17g must match main_sim time %.17g at registration",
                                   sub_sim->t, bridge->main_sim->t);
    }
    if (validate_subsystem_canonical_state(&bridge->main_sim->particles[host_index], sub_sim,
                                          check_barycenter) != 0) {
        return -1;
    }
    if (!isfinite(dt_inner) || dt_inner <= 0.0) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "dt_inner must be finite and positive");
    if (dt_inner > bridge->dt_outer) return bridge_error_code(REB_BRIDGE_ERROR_INVALID_ARGUMENT, "dt_inner must be <= dt_outer");

    double ratio = bridge->dt_outer / dt_inner;
    if (!isfinite(ratio) || fabs(ratio - round(ratio)) > 1.0e-12) {
        return bridge_error_format(REB_BRIDGE_ERROR_INVALID_ARGUMENT,
                                   "dt_outer/dt_inner must be a finite integer ratio (received %.17g/%.17g = %.17g)",
                                   bridge->dt_outer, dt_inner, ratio);
    }

    if (bridge->n_subsystems == bridge->cap_subsystems) {
        int new_cap = bridge->cap_subsystems == 0 ? 4 : bridge->cap_subsystems * 2;
        struct reb_bridge_subsystem* new_subsystems =
            (struct reb_bridge_subsystem*)realloc(
                bridge->subsystems,
                sizeof(*new_subsystems) * (size_t)new_cap
            );
        if (!new_subsystems) return bridge_error_code(REB_BRIDGE_ERROR_OUT_OF_MEMORY, "cannot grow subsystem registry");
        bridge->subsystems = new_subsystems;
        bridge->cap_subsystems = new_cap;
    }

    sub_sim->dt = dt_inner;
    struct reb_bridge_subsystem* sub = &bridge->subsystems[bridge->n_subsystems++];
    memset(sub, 0, sizeof(*sub));
    sub->sub_sim = sub_sim;
    sub->host_index = host_index;
    sub->dt_inner = dt_inner;
    sub->owns_sub_sim = owns_sub_sim ? 1 : 0;
    return 0;
}

int reb_bridge_add_subsystem(
    struct reb_bridge* bridge,
    int host_index,
    struct reb_simulation* sub_sim,
    double dt_inner,
    int owns_sub_sim
) {
    return reb_bridge_internal_add_subsystem(bridge, host_index, sub_sim, dt_inner,
                                            owns_sub_sim, 1);
}

int reb_bridge_subsystem_count(const struct reb_bridge* bridge) {
    return bridge ? bridge->n_subsystems : 0;
}

struct reb_simulation* reb_bridge_main_sim(struct reb_bridge* bridge) {
    return bridge ? bridge->main_sim : NULL;
}

struct reb_simulation* reb_bridge_sub_sim(struct reb_bridge* bridge, int subsystem_index) {
    if (!bridge || subsystem_index < 0 || subsystem_index >= bridge->n_subsystems) return NULL;
    return bridge->subsystems[subsystem_index].sub_sim;
}
