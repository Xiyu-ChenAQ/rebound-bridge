/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef VAL_PARTICLE_H
#define VAL_PARTICLE_H

/**
 * @file  val_particle.h
 * @brief Inline particle geometry helpers shared across validation programs.
 */

#include <math.h>
#include <string.h>
#include "rebound.h"

static inline struct reb_particle val_particle_from_state(
    double m,
    double x,  double y,  double z,
    double vx, double vy, double vz
) {
    struct reb_particle p;
    memset(&p, 0, sizeof(p));
    p.m = m;
    p.x = x;   p.y = y;   p.z = z;
    p.vx = vx; p.vy = vy; p.vz = vz;
    return p;
}

/* Return a copy of rel with host's position/velocity added. */
static inline struct reb_particle val_particle_add(
    const struct reb_particle* host,
    const struct reb_particle* rel
) {
    struct reb_particle p = *rel;
    p.x  += host->x;  p.y  += host->y;  p.z  += host->z;
    p.vx += host->vx; p.vy += host->vy; p.vz += host->vz;
    return p;
}

static inline double val_distance(
    const struct reb_particle* a,
    const struct reb_particle* b
) {
    const double dx = a->x - b->x;
    const double dy = a->y - b->y;
    const double dz = a->z - b->z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/*
 * ||(a_sec - a_pri) - (b_sec - b_pri)||  - relative position error between
 * two integrators for the same pair.
 */
static inline double val_relative_position_delta(
    const struct reb_particle* a_pri, const struct reb_particle* a_sec,
    const struct reb_particle* b_pri, const struct reb_particle* b_sec
) {
    const double dx = (a_sec->x - a_pri->x) - (b_sec->x - b_pri->x);
    const double dy = (a_sec->y - a_pri->y) - (b_sec->y - b_pri->y);
    const double dz = (a_sec->z - a_pri->z) - (b_sec->z - b_pri->z);
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/*
 * Mass-weighted barycenter of particles at the given indices array.
 * Returns a zeroed particle if total mass is zero.
 */
static inline struct reb_particle val_barycenter(
    const struct reb_particle* bodies,
    const int* indices,
    int count
) {
    struct reb_particle com;
    memset(&com, 0, sizeof(com));
    for (int i = 0; i < count; i++) {
        const struct reb_particle* p = &bodies[indices[i]];
        com.m  += p->m;
        com.x  += p->m * p->x;  com.y  += p->m * p->y;  com.z  += p->m * p->z;
        com.vx += p->m * p->vx; com.vy += p->m * p->vy; com.vz += p->m * p->vz;
    }
    if (com.m > 0.0) {
        com.x  /= com.m; com.y  /= com.m; com.z  /= com.m;
        com.vx /= com.m; com.vy /= com.m; com.vz /= com.m;
    }
    return com;
}

#endif /* VAL_PARTICLE_H */
