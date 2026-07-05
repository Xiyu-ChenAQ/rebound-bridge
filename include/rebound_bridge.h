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
 *   EOS_LF4_2   - 4th-order EOS-LF composition (Blanes & Casas).
 *
 * Experimental / validation targets only - implementations live in
 * src/bridge_fg4_experimental.c and src/bridge_force_gradient_experimental.c:
 *   EOS_LF8_6_4 - 8th-order EOS-LF composition; extra substeps, not faster in practice.
 *   FORWARD_FG4 - 4th-order force-gradient forward scheme; requires second
 *                 force-gradient evaluation per step.
 */
enum reb_bridge_method {
    REB_BRIDGE_METHOD_INVALID     = -1,
    REB_BRIDGE_METHOD_ORDER2      =  0,
    REB_BRIDGE_METHOD_EOS_LF4_2   =  1,
    REB_BRIDGE_METHOD_EOS_LF8_6_4 =  2, /* experimental */
    REB_BRIDGE_METHOD_FORWARD_FG4 =  3  /* experimental */
};

enum reb_bridge_cross_kick_kernel {
    REB_BRIDGE_CROSS_KICK_INVALID        = -1,
    REB_BRIDGE_CROSS_KICK_PLAIN          =  0,
    REB_BRIDGE_CROSS_KICK_FORCE_GRADIENT =  1  /* experimental */
};

enum reb_bridge_corrector {
    REB_BRIDGE_CORRECTOR_INVALID       = -1,
    REB_BRIDGE_CORRECTOR_NONE          =  0,
    REB_BRIDGE_CORRECTOR_COMMUTATOR_AB =  1  /* experimental */
};

enum reb_bridge_corrector_direction {
    REB_BRIDGE_CORRECTOR_INVERSE = -1,
    REB_BRIDGE_CORRECTOR_FORWARD =  1
};

struct reb_bridge_force_gradient_check {
    double max_abs_difference;
    double max_relative_difference;
    double max_abs_analytic;
    double max_abs_finite_difference;
    double internal_weighted_acc_residual;
    double main_weighted_acc_residual;
    double force_gradient_internal_weighted_acc_residual;
    double force_gradient_main_weighted_acc_residual;
};


REB_BRIDGE_API struct reb_bridge* reb_bridge_create(struct reb_simulation* main_sim, double dt_outer);
REB_BRIDGE_API void reb_bridge_set_owns_main(struct reb_bridge* bridge, int owns_main_sim);
REB_BRIDGE_API int reb_bridge_set_method(struct reb_bridge* bridge, enum reb_bridge_method method);
REB_BRIDGE_API enum reb_bridge_method reb_bridge_get_method(const struct reb_bridge* bridge);
REB_BRIDGE_API void reb_bridge_free(struct reb_bridge* bridge);

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
REB_BRIDGE_API int reb_bridge_step(struct reb_bridge* bridge, double dt);
REB_BRIDGE_API int reb_bridge_advance(struct reb_bridge* bridge, double duration);
REB_BRIDGE_API int reb_bridge_advance_eos_lf4_2(struct reb_bridge* bridge, double duration);

REB_BRIDGE_API int reb_bridge_set_cross_kick_kernel(
    struct reb_bridge* bridge,
    enum reb_bridge_cross_kick_kernel kernel
);
REB_BRIDGE_API enum reb_bridge_cross_kick_kernel reb_bridge_get_cross_kick_kernel(const struct reb_bridge* bridge);
REB_BRIDGE_API int reb_bridge_set_force_gradient_scale(struct reb_bridge* bridge, double scale);
REB_BRIDGE_API double reb_bridge_get_force_gradient_scale(const struct reb_bridge* bridge);

REB_BRIDGE_API int reb_bridge_set_corrector(struct reb_bridge* bridge, enum reb_bridge_corrector corrector);
REB_BRIDGE_API enum reb_bridge_corrector reb_bridge_get_corrector(const struct reb_bridge* bridge);
REB_BRIDGE_API int reb_bridge_set_corrector_scale(struct reb_bridge* bridge, double scale);
REB_BRIDGE_API double reb_bridge_get_corrector_scale(const struct reb_bridge* bridge);
REB_BRIDGE_API int reb_bridge_apply_corrector(
    struct reb_bridge* bridge,
    enum reb_bridge_corrector_direction direction
);

REB_BRIDGE_API int reb_bridge_force_gradient_finite_difference_check(
    struct reb_bridge* bridge,
    int subsystem_index,
    double epsilon,
    struct reb_bridge_force_gradient_check* out
);

REB_BRIDGE_API int reb_bridge_set_integrator_whfast(struct reb_simulation* sim);
REB_BRIDGE_API int reb_bridge_set_integrator_ias15(struct reb_simulation* sim);

#ifdef __cplusplus
}
#endif

#endif /* REBOUND_BRIDGE_H */
