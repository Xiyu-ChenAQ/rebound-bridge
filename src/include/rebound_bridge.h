/* SPDX-License-Identifier: GPL-3.0-or-later */

#ifndef REBOUND_BRIDGE_H
#define REBOUND_BRIDGE_H

/**
 * @file    rebound_bridge.h
 * @brief   Standalone symplectic bridge scheduler for REBOUND simulations.
 */

#include "rebound.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef REB_BRIDGE_API
#if defined(_WIN32) && defined(REB_BRIDGE_BUILD_SHARED)
#if defined(REB_BRIDGE_BUILDING_LIBRARY)
#define REB_BRIDGE_API __declspec(dllexport)
#else
#define REB_BRIDGE_API __declspec(dllimport)
#endif
#else
#define REB_BRIDGE_API
#endif
#endif

#define REB_BRIDGE_G_AU_YR_MSUN (4.0 * M_PI * M_PI)

struct reb_bridge;

/**
 * Bridge operator compositions.
 *
 * Stable for production use:
 *   ORDER2      - 2nd-order symmetric (Verlet-style) bridge.
 *   EOS_LF4_2   - symmetric EOS-LF composition of generalized order (4,2):
 *                 for H = A + epsilon B, its global splitting error over a
 *                 fixed integration interval has leading structure
 *                 O(epsilon h^4 + epsilon^2 h^2).
 *
 */
enum reb_bridge_method {
    REB_BRIDGE_METHOD_INVALID     = -1,
    REB_BRIDGE_METHOD_ORDER2      =  0,
    REB_BRIDGE_METHOD_EOS_LF4_2   =  1,
};

enum reb_bridge_subsystem_coupling {
    REB_BRIDGE_SUBSYSTEM_COUPLING_INVALID  = -1,
    REB_BRIDGE_SUBSYSTEM_COUPLING_MONOPOLE =  0,
    REB_BRIDGE_SUBSYSTEM_COUPLING_RESOLVED =  1
};

/** Create a bridge with a finite positive outer timestep. The main simulation
 * must have finite time, a finite non-negative G, and finite particle masses,
 * positions and velocities. Particle masses must be non-negative. */
REB_BRIDGE_API struct reb_bridge* reb_bridge_create(struct reb_simulation* main_sim, double dt_outer);
REB_BRIDGE_API void reb_bridge_set_owns_main(struct reb_bridge* bridge, int owns_main_sim);
REB_BRIDGE_API int reb_bridge_set_method(struct reb_bridge* bridge, enum reb_bridge_method method);
REB_BRIDGE_API enum reb_bridge_method reb_bridge_get_method(const struct reb_bridge* bridge);
REB_BRIDGE_API void reb_bridge_free(struct reb_bridge* bridge);

/**
 * Register a resolved subsystem represented by a host particle in the main
 * simulation.  Individual subsystem particle masses must be finite and
 * non-negative.  The subsystem total mass and matching host mass must be
 * finite, strictly positive, and equal.  Subsystem coordinates and velocities
 * must be barycentric (zero center-of-mass position and velocity).
 * Main and subsystem simulations must use the same G and initial time, and
 * must be distinct objects. dt_inner must be finite and positive, with an
 * integer dt_outer/dt_inner ratio. Registration does not transfer ownership
 * on failure.
 * For fixed-step integrators, dt_inner is the nominal maximum substep; each
 * operator-composition segment is divided into equal steps no larger than it.
 */
REB_BRIDGE_API int reb_bridge_add_subsystem(
    struct reb_bridge* bridge,
    int host_index,
    struct reb_simulation* sub_sim,
    double dt_inner,
    int owns_sub_sim
);

REB_BRIDGE_API int reb_bridge_subsystem_count(const struct reb_bridge* bridge);
REB_BRIDGE_API struct reb_simulation* reb_bridge_main_sim(struct reb_bridge* bridge);
REB_BRIDGE_API struct reb_simulation* reb_bridge_sub_sim(struct reb_bridge* bridge, int subsystem_index);

REB_BRIDGE_API int reb_bridge_apply_cross_kick(struct reb_bridge* bridge, double dt_kick);
/**
 * Advance one Bridge step of size dt using the method selected with
 * reb_bridge_set_method().  New bridges default to ORDER2.
 */
REB_BRIDGE_API int reb_bridge_step(struct reb_bridge* bridge, double dt);
/** Advance by a finite non-negative duration containing an integer number of
 * outer timesteps. A zero duration is a valid no-op. */
REB_BRIDGE_API int reb_bridge_advance(struct reb_bridge* bridge, double duration);
REB_BRIDGE_API int reb_bridge_advance_eos_lf4_2(struct reb_bridge* bridge, double duration);
REB_BRIDGE_API int reb_bridge_set_subsystem_coupling(
    struct reb_bridge* bridge,
    enum reb_bridge_subsystem_coupling coupling
);
REB_BRIDGE_API enum reb_bridge_subsystem_coupling reb_bridge_get_subsystem_coupling(
    const struct reb_bridge* bridge
);

/** Enable (1, the default) or disable (0) WHFast Jacobi-state reuse and direct
 * Jacobi velocity updates. Disabling retains Cartesian residual increments
 * and the component's configured safe_mode. Switching modes invalidates all
 * component caches without changing the current physical state. */
REB_BRIDGE_API int reb_bridge_set_whfast_jacobi_updates(struct reb_bridge* bridge, int enabled);
/** Return 0 or 1, or -1 for NULL. Like other getters, preserves error state. */
REB_BRIDGE_API int reb_bridge_get_whfast_jacobi_updates(const struct reb_bridge* bridge);

REB_BRIDGE_API int reb_bridge_set_integrator_whfast(struct reb_simulation* sim);
REB_BRIDGE_API int reb_bridge_set_integrator_ias15(struct reb_simulation* sim);

enum reb_bridge_error {
    REB_BRIDGE_ERROR_NONE = 0,
    REB_BRIDGE_ERROR_INVALID_ARGUMENT = 1,
    REB_BRIDGE_ERROR_OUT_OF_MEMORY = 2,
    REB_BRIDGE_ERROR_DURATION_NOT_DIVISIBLE = 3,
    REB_BRIDGE_ERROR_DUPLICATE_REGISTRATION = 4,
    REB_BRIDGE_ERROR_INTERNAL = 5
};

/** Error state is local to the calling thread. Mutating public operations clear
 * it before attempting work; getters and reb_bridge_free() preserve it.
 * The message pointer remains valid until the next mutation or explicit clear
 * on this thread. An empty string accompanies REB_BRIDGE_ERROR_NONE. */
REB_BRIDGE_API enum reb_bridge_error reb_bridge_last_error(void);
REB_BRIDGE_API const char* reb_bridge_last_error_message(void);
REB_BRIDGE_API void reb_bridge_clear_error(void);

#ifdef __cplusplus
}
#endif

#endif /* REBOUND_BRIDGE_H */
