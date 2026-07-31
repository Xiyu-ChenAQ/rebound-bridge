/* SPDX-License-Identifier: GPL-3.0-or-later */
/* Xiyu-Chen */
#ifndef REBOUND_BRIDGE_INTERNAL_H
#define REBOUND_BRIDGE_INTERNAL_H

#include "rebound_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct reb_bridge_vec3 {
    double x;
    double y;
    double z;
};

struct reb_bridge_subsystem {
    struct reb_simulation* sub_sim;
    int host_index;
    double dt_inner;
    int owns_sub_sim;
};

struct reb_bridge {
    struct reb_simulation* main_sim;
    double dt_outer;
    int owns_main_sim;
    enum reb_bridge_method method;
    enum reb_bridge_cross_kick_kernel cross_kick_kernel;
    enum reb_bridge_subsystem_coupling subsystem_coupling;
    double force_gradient_scale;
    enum reb_bridge_corrector corrector;
    double corrector_scale;
    struct reb_bridge_subsystem* subsystems;
    int n_subsystems;
    int cap_subsystems;
    struct reb_bridge_vec3* internal_acc;
    struct reb_bridge_vec3* source_acc;
    struct reb_bridge_vec3* main_acc;
    struct reb_bridge_vec3* fg_internal_acc;
    struct reb_bridge_vec3* fg_main_acc;
    struct reb_bridge_vec3* fd_internal_acc;
    struct reb_bridge_vec3* fd_main_acc;
    struct reb_bridge_vec3* saved_sub_pos;
    struct reb_bridge_vec3* saved_main_pos;
    int cap_sub_particles;
    int cap_main_particles;
};

static inline struct reb_bridge_vec3 vec3(double x, double y, double z) {
    struct reb_bridge_vec3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

static inline struct reb_bridge_vec3 particle_pos(const struct reb_particle* p) {
    return vec3(p->x, p->y, p->z);
}

static inline struct reb_bridge_vec3 vec_add(struct reb_bridge_vec3 a, struct reb_bridge_vec3 b) {
    return vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline struct reb_bridge_vec3 vec_sub(struct reb_bridge_vec3 a, struct reb_bridge_vec3 b) {
    return vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline struct reb_bridge_vec3 vec_scale(struct reb_bridge_vec3 v, double s) {
    return vec3(v.x * s, v.y * s, v.z * s);
}

static inline double vec_dot(struct reb_bridge_vec3 a, struct reb_bridge_vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline double vec_norm(struct reb_bridge_vec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

static inline double bridge_time_sync_tolerance(double a, double b) {
    const double scale = fmax(fabs(a), fabs(b));
    return 1.0e-12 + 1.0e-15 * fmax(1.0, scale);
}

static inline double bridge_time_advance_tolerance(double a, double b) {
    const double scale = fmax(fabs(a), fabs(b));
    return 1.0e-14 * fmax(1.0, scale);
}

static inline int bridge_reached_time(double t, double target) {
    return t >= target - bridge_time_advance_tolerance(t, target);
}

static inline int bridge_error(const char* msg) {
    fprintf(stderr, "[rebound_bridge] %s\n", msg);
    return -1;
}

static inline int valid_particle_index(const struct reb_simulation* sim, int index) {
    return sim && index >= 0 && index < sim->N;
}

static inline int subsystem_barycenter(
    const struct reb_simulation* sim,
    struct reb_bridge_vec3* out
) {
    if (sim->N <= 0) return bridge_error("subsystem has no particles");

    double mt = 0.0;
    struct reb_bridge_vec3 weighted = vec3(0.0, 0.0, 0.0);
    for (int i = 0; i < sim->N; i++) {
        const struct reb_particle* p = &sim->particles[i];
        mt += p->m;
        weighted = vec_add(weighted, vec_scale(particle_pos(p), p->m));
    }
    if (mt <= 0.0) return bridge_error("subsystem total mass must be positive");

    *out = vec_scale(weighted, 1.0 / mt);
    return 0;
}

static inline struct reb_bridge_vec3 gravity_from_source(
    struct reb_bridge_vec3 pos,
    const struct reb_particle* source,
    double G
) {
    struct reb_bridge_vec3 r = vec_sub(pos, particle_pos(source));
    double rnorm = vec_norm(r);
    double inv_r3 = 1.0 / (rnorm * rnorm * rnorm);
    return vec_scale(r, -G * source->m * inv_r3);
}

static inline struct reb_bridge_vec3 gravity_from_point_mass(
    struct reb_bridge_vec3 pos,
    struct reb_bridge_vec3 source_pos,
    double source_mass,
    double G
) {
    struct reb_bridge_vec3 r = vec_sub(pos, source_pos);
    double rnorm = vec_norm(r);
    double inv_r3 = 1.0 / (rnorm * rnorm * rnorm);
    return vec_scale(r, -G * source_mass * inv_r3);
}

int reb_bridge_internal_ensure_kick_workspace(struct reb_bridge* bridge, int sub_particles, int main_particles);
int reb_bridge_internal_check_subsystem_time_sync(struct reb_bridge* bridge, struct reb_bridge_subsystem* sub);
int reb_bridge_internal_compute_subsystem_cross_acc(
    struct reb_bridge* bridge,
    struct reb_bridge_subsystem* sub,
    struct reb_bridge_vec3* internal_acc,
    struct reb_bridge_vec3* main_acc
);
int reb_bridge_internal_apply_force_gradient_correction(
    struct reb_bridge* bridge,
    struct reb_bridge_subsystem* sub,
    struct reb_bridge_vec3* internal_acc,
    struct reb_bridge_vec3* main_acc,
    double step_dt
);
int reb_bridge_internal_apply_cross_kick_for_step(struct reb_bridge* bridge, double dt_kick, double step_dt);
int reb_bridge_internal_apply_cross_kick_with_kernel_for_step(
    struct reb_bridge* bridge,
    double dt_kick,
    double step_dt,
    enum reb_bridge_cross_kick_kernel kernel
);
int reb_bridge_internal_integrate_all_to_target(struct reb_bridge* bridge, double target);
int reb_bridge_internal_integrate_all_to_target_symmetric(struct reb_bridge* bridge, double target);
void reb_bridge_internal_snap_all_times(struct reb_bridge* bridge, double time);
int reb_bridge_internal_step_order2(struct reb_bridge* bridge, double dt);
int reb_bridge_internal_step_eos_lf4_2(struct reb_bridge* bridge, double dt);
int reb_bridge_internal_step_eos_lf8_6_4(struct reb_bridge* bridge, double dt);
int reb_bridge_internal_step_forward_fg4(struct reb_bridge* bridge, double dt);

/* Experimental advance functions - implemented in bridge_fg4_experimental.c. */
int reb_bridge_advance_eos_lf8_6_4(struct reb_bridge* bridge, double duration);
int reb_bridge_advance_forward_fg4(struct reb_bridge* bridge, double duration);

#endif /* REBOUND_BRIDGE_INTERNAL_H */
