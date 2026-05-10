#ifndef REBOUND_BRIDGE_H
#define REBOUND_BRIDGE_H

/**
 * @file    rebound_bridge.h
 * @brief   Standalone symplectic bridge scheduler for REBOUND simulations.
 */

#include "rebound.h"

#ifdef __cplusplus
extern "C" {
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

DLLEXPORT struct reb_bridge* reb_bridge_create(struct reb_simulation* main_sim, double dt_outer);
DLLEXPORT void reb_bridge_set_owns_main(struct reb_bridge* bridge, int owns_main_sim);
DLLEXPORT void reb_bridge_free(struct reb_bridge* bridge);

DLLEXPORT int reb_bridge_add_subsystem(
    struct reb_bridge* bridge,
    int host_index,
    struct reb_simulation* sub_sim,
    double dt_inner,
    int owns_sub_sim
);

DLLEXPORT int reb_bridge_subsystem_count(const struct reb_bridge* bridge);
DLLEXPORT struct reb_simulation* reb_bridge_main_sim(struct reb_bridge* bridge);
DLLEXPORT struct reb_simulation* reb_bridge_sub_sim(struct reb_bridge* bridge, int subsystem_index);

DLLEXPORT int reb_bridge_apply_cross_kick(struct reb_bridge* bridge, double dt_kick);
DLLEXPORT int reb_bridge_step(struct reb_bridge* bridge, double dt);
DLLEXPORT int reb_bridge_advance(struct reb_bridge* bridge, double duration);
DLLEXPORT int reb_bridge_snapshot(
    const struct reb_bridge* bridge,
    int subsystem_index,
    struct reb_bridge_snapshot* out
);

DLLEXPORT struct reb_simulation* reb_bridge_make_main_sun_emb(void);
DLLEXPORT struct reb_simulation* reb_bridge_make_sub_earth_moon(void);
DLLEXPORT struct reb_bridge* reb_bridge_create_earth_moon(double dt_outer, int sub_ratio);

#ifdef __cplusplus
}
#endif

#endif /* REBOUND_BRIDGE_H */
