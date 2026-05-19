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

struct reb_bridge_body_state {
    double x, y, z;
    double vx, vy, vz;
    double m;
};

struct reb_bridge_snapshot {
    double time;
    int host_index;
    int member_count;
    double host_primary_distance;
    double subsystem_barycenter_distance;
    double subsystem_internal_distance;
    struct reb_bridge_body_state host;
    struct reb_bridge_body_state member0;
    struct reb_bridge_body_state member1;
};

REB_BRIDGE_API struct reb_bridge* reb_bridge_create(struct reb_simulation* main_sim, double dt_outer);
REB_BRIDGE_API void reb_bridge_set_owns_main(struct reb_bridge* bridge, int owns_main_sim);
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
REB_BRIDGE_API int reb_bridge_snapshot(
    const struct reb_bridge* bridge,
    int subsystem_index,
    struct reb_bridge_snapshot* out
);

REB_BRIDGE_API int reb_bridge_set_integrator_whfast(struct reb_simulation* sim);
REB_BRIDGE_API int reb_bridge_set_integrator_ias15(struct reb_simulation* sim);

REB_BRIDGE_API struct reb_simulation* reb_bridge_make_main_sun_emb(void);
REB_BRIDGE_API struct reb_simulation* reb_bridge_make_sub_earth_moon(void);
REB_BRIDGE_API struct reb_bridge* reb_bridge_create_earth_moon(double dt_outer, int sub_ratio);

#ifdef __cplusplus
}
#endif

#endif /* REBOUND_BRIDGE_H */
