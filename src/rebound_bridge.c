/* SPDX-License-Identifier: GPL-3.0-only */

#include "rebound_bridge.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    struct reb_bridge_subsystem* subsystems;
    int n_subsystems;
    int cap_subsystems;
    struct reb_vec3d* internal_acc;
    struct reb_vec3d* source_acc;
    struct reb_vec3d* main_acc;
    int cap_sub_particles;
    int cap_main_particles;
};

static double bridge_time_sync_tolerance(double a, double b) {
    const double scale = fmax(fabs(a), fabs(b));
    return 1.0e-12 + 1.0e-15 * fmax(1.0, scale);
}

static int bridge_error(const char* msg) {
    fprintf(stderr, "[rebound_bridge] %s\n", msg);
    return -1;
}

int reb_bridge_set_integrator_whfast(struct reb_simulation* sim) {
    struct reb_integrator_whfast_state* whfast =
        (struct reb_integrator_whfast_state*)reb_simulation_set_integrator(sim, "whfast");
    if (!whfast) return bridge_error("failed to set REBOUND integrator to whfast");
    whfast->safe_mode = 1;
    return 0;
}

int reb_bridge_set_integrator_saba(struct reb_simulation* sim, int type) {
    struct reb_integrator_saba_state* saba =
        (struct reb_integrator_saba_state*)reb_simulation_set_integrator(sim, "saba");
    if (!saba) return bridge_error("failed to set REBOUND integrator to saba");
    if (type != 0) saba->type = type;
    saba->safe_mode = 1;
    return 0;
}

static int valid_particle_index(const struct reb_simulation* sim, int index) {
    return sim && index >= 0 && (size_t)index < sim->N;
}

static struct reb_bridge_body_state body_state(const struct reb_particle* p) {
    struct reb_bridge_body_state s;
    s.x = p->x;
    s.y = p->y;
    s.z = p->z;
    s.vx = p->vx;
    s.vy = p->vy;
    s.vz = p->vz;
    s.m = p->m;
    return s;
}

static double particle_distance(const struct reb_particle* a, const struct reb_particle* b) {
    const double dx = a->x - b->x;
    const double dy = a->y - b->y;
    const double dz = a->z - b->z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

static int subsystem_barycenter(
    const struct reb_simulation* sim,
    struct reb_vec3d* out
) {
    if (sim->N <= 0) return bridge_error("subsystem has no particles");

    double mt = 0.0;
    double x = 0.0, y = 0.0, z = 0.0;
    for (size_t i = 0; i < sim->N; i++) {
        const struct reb_particle* p = &sim->particles[i];
        mt += p->m;
        x += p->x * p->m;
        y += p->y * p->m;
        z += p->z * p->m;
    }
    if (mt <= 0.0) return bridge_error("subsystem total mass must be positive");

    out->x = x / mt;
    out->y = y / mt;
    out->z = z / mt;
    return 0;
}

/* Gravitational acceleration at (px,py,pz) due to point mass m at (sx,sy,sz). */
static struct reb_vec3d accel_from_mass(
    double px, double py, double pz,
    double sx, double sy, double sz,
    double m, double G
) {
    const double rx = px - sx;
    const double ry = py - sy;
    const double rz = pz - sz;
    const double r = sqrt(rx * rx + ry * ry + rz * rz);
    const double inv_r3 = 1.0 / (r * r * r);
    const double prefact = -G * m * inv_r3;
    struct reb_vec3d a = { rx * prefact, ry * prefact, rz * prefact };
    return a;
}

struct reb_bridge* reb_bridge_create(struct reb_simulation* main_sim, double dt_outer) {
    if (dt_outer <= 0.0) {
        bridge_error("dt_outer must be positive");
        return NULL;
    }

    struct reb_bridge* bridge = (struct reb_bridge*)calloc(1, sizeof(*bridge));
    if (!bridge) {
        bridge_error("out of memory");
        return NULL;
    }

    bridge->main_sim = main_sim;
    bridge->dt_outer = dt_outer;
    bridge->owns_main_sim = 0;
    bridge->main_sim->dt = dt_outer;
    return bridge;
}

void reb_bridge_set_owns_main(struct reb_bridge* bridge, int owns_main_sim) {
    bridge->owns_main_sim = owns_main_sim;
}

void reb_bridge_free(struct reb_bridge* bridge) {
    if (!bridge) return;
    for (int i = 0; i < bridge->n_subsystems; i++) {
        if (bridge->subsystems[i].owns_sub_sim && bridge->subsystems[i].sub_sim) {
            reb_simulation_free(bridge->subsystems[i].sub_sim);
        }
    }
    free(bridge->subsystems);
    free(bridge->internal_acc);
    free(bridge->source_acc);
    free(bridge->main_acc);
    if (bridge->owns_main_sim && bridge->main_sim) {
        reb_simulation_free(bridge->main_sim);
    }
    free(bridge);
}

int reb_bridge_add_subsystem(
    struct reb_bridge* bridge,
    int host_index,
    struct reb_simulation* sub_sim,
    double dt_inner,
    int owns_sub_sim
) {
    if (!valid_particle_index(bridge->main_sim, host_index)) {
        return bridge_error("host_index is out of bounds");
    }
    if (sub_sim->N <= 0) return bridge_error("sub_sim has no particles");
    if (dt_inner <= 0.0) return bridge_error("dt_inner must be positive");
    if (dt_inner > bridge->dt_outer) return bridge_error("dt_inner must be <= dt_outer");

    double ratio = bridge->dt_outer / dt_inner;
    if (fabs(ratio - round(ratio)) > 1.0e-12) {
        return bridge_error("dt_outer/dt_inner must be an integer ratio");
    }

    if (bridge->n_subsystems == bridge->cap_subsystems) {
        int new_cap = bridge->cap_subsystems == 0 ? 4 : bridge->cap_subsystems * 2;
        struct reb_bridge_subsystem* new_subsystems =
            (struct reb_bridge_subsystem*)realloc(
                bridge->subsystems,
                sizeof(*new_subsystems) * (size_t)new_cap
            );
        if (!new_subsystems) return bridge_error("out of memory");
        bridge->subsystems = new_subsystems;
        bridge->cap_subsystems = new_cap;
    }

    sub_sim->dt = dt_inner;
    struct reb_bridge_subsystem* sub = &bridge->subsystems[bridge->n_subsystems++];
    sub->sub_sim = sub_sim;
    sub->host_index = host_index;
    sub->dt_inner = dt_inner;
    sub->owns_sub_sim = owns_sub_sim;
    return 0;
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

static int ensure_kick_workspace(struct reb_bridge* bridge, int sub_particles, int main_particles) {
    if (sub_particles > bridge->cap_sub_particles) {
        struct reb_vec3d* internal_acc =
            (struct reb_vec3d*)realloc(
                bridge->internal_acc,
                sizeof(*internal_acc) * (size_t)sub_particles
            );
        if (!internal_acc) return bridge_error("out of memory");
        bridge->internal_acc = internal_acc;

        struct reb_vec3d* source_acc =
            (struct reb_vec3d*)realloc(
                bridge->source_acc,
                sizeof(*source_acc) * (size_t)sub_particles
            );
        if (!source_acc) return bridge_error("out of memory");
        bridge->source_acc = source_acc;
        bridge->cap_sub_particles = sub_particles;
    }

    if (main_particles > bridge->cap_main_particles) {
        struct reb_vec3d* main_acc =
            (struct reb_vec3d*)realloc(
                bridge->main_acc,
                sizeof(*main_acc) * (size_t)main_particles
            );
        if (!main_acc) return bridge_error("out of memory");
        bridge->main_acc = main_acc;
        bridge->cap_main_particles = main_particles;
    }

    return 0;
}

static int apply_subsystem_cross_kick(
    struct reb_bridge* bridge,
    struct reb_bridge_subsystem* sub,
    double dt_kick
) {
    struct reb_simulation* main_sim = bridge->main_sim;
    struct reb_simulation* sub_sim = sub->sub_sim;
    if (!valid_particle_index(main_sim, sub->host_index)) {
        return bridge_error("host_index is out of bounds during kick");
    }
    const double time_diff = sub_sim->t - main_sim->t;
    const double time_tol = bridge_time_sync_tolerance(main_sim->t, sub_sim->t);
    if (fabs(time_diff) > time_tol) {
        return bridge_error("main and subsystem times are out of sync");
    }
    /* Snap only the time label inside floating-point tolerance; particle state is unchanged. */
    sub_sim->t = main_sim->t;

    const int n = sub_sim->N;
    struct reb_vec3d bc;
    if (subsystem_barycenter(sub_sim, &bc) != 0) return -1;

    if (ensure_kick_workspace(bridge, n, main_sim->N) != 0) return -1;
    struct reb_vec3d* internal_acc = bridge->internal_acc;
    struct reb_vec3d* source_acc = bridge->source_acc;
    struct reb_vec3d* main_acc = bridge->main_acc;
    memset(internal_acc, 0, sizeof(*internal_acc) * (size_t)n);
    memset(main_acc, 0, sizeof(*main_acc) * (size_t)main_sim->N);

    const struct reb_particle* host = &main_sim->particles[sub->host_index];
    const double bc_x = host->x + bc.x;
    const double bc_y = host->y + bc.y;
    const double bc_z = host->z + bc.z;
    const double G = main_sim->G;
    double mt = 0.0;
    for (int i = 0; i < n; i++) {
        mt += sub_sim->particles[i].m;
    }
    if (mt <= 0.0) {
        return bridge_error("subsystem total mass must be positive");
    }

    double host_ax = 0.0, host_ay = 0.0, host_az = 0.0;

    /* Apply the external tidal field to subsystem members; feed the equal reaction back to the main system. */
    for (int source_index = 0; source_index < (int)main_sim->N; source_index++) {
        if (source_index == sub->host_index) continue;
        const struct reb_particle* source = &main_sim->particles[source_index];
        if (source->m == 0.0) continue;

        struct reb_vec3d a_mono = accel_from_mass(bc_x, bc_y, bc_z, source->x, source->y, source->z, source->m, G);
        double com_x = 0.0, com_y = 0.0, com_z = 0.0;
        double src_x = 0.0, src_y = 0.0, src_z = 0.0;
        for (int i = 0; i < n; i++) {
            const struct reb_particle* member = &sub_sim->particles[i];
            const double gx = host->x + member->x;
            const double gy = host->y + member->y;
            const double gz = host->z + member->z;
            source_acc[i] = accel_from_mass(gx, gy, gz, source->x, source->y, source->z, source->m, G);
            com_x += source_acc[i].x * member->m;
            com_y += source_acc[i].y * member->m;
            com_z += source_acc[i].z * member->m;
            struct reb_vec3d a_src = accel_from_mass(source->x, source->y, source->z, gx, gy, gz, member->m, G);
            src_x += a_src.x;
            src_y += a_src.y;
            src_z += a_src.z;
        }

        const double inv_mt = 1.0 / mt;
        com_x *= inv_mt;
        com_y *= inv_mt;
        com_z *= inv_mt;
        host_ax += com_x - a_mono.x;
        host_ay += com_y - a_mono.y;
        host_az += com_z - a_mono.z;

        for (int i = 0; i < n; i++) {
            internal_acc[i].x += source_acc[i].x - com_x;
            internal_acc[i].y += source_acc[i].y - com_y;
            internal_acc[i].z += source_acc[i].z - com_z;
        }

        struct reb_vec3d a_src_mono = accel_from_mass(source->x, source->y, source->z, bc_x, bc_y, bc_z, mt, G);
        main_acc[source_index].x += src_x - a_src_mono.x;
        main_acc[source_index].y += src_y - a_src_mono.y;
        main_acc[source_index].z += src_z - a_src_mono.z;
    }

    for (int i = 0; i < n; i++) {
        struct reb_particle* p = &sub_sim->particles[i];
        p->vx += internal_acc[i].x * dt_kick;
        p->vy += internal_acc[i].y * dt_kick;
        p->vz += internal_acc[i].z * dt_kick;
    }

    struct reb_particle* host_mut = &main_sim->particles[sub->host_index];
    host_mut->vx += host_ax * dt_kick;
    host_mut->vy += host_ay * dt_kick;
    host_mut->vz += host_az * dt_kick;

    for (int source_index = 0; source_index < (int)main_sim->N; source_index++) {
        if (source_index == sub->host_index) continue;
        struct reb_particle* source = &main_sim->particles[source_index];
        source->vx += main_acc[source_index].x * dt_kick;
        source->vy += main_acc[source_index].y * dt_kick;
        source->vz += main_acc[source_index].z * dt_kick;
    }

    reb_simulation_synchronize(sub_sim);
    reb_simulation_synchronize(main_sim);
    return 0;
}

int reb_bridge_apply_cross_kick(struct reb_bridge* bridge, double dt_kick) {
    for (int i = 0; i < bridge->n_subsystems; i++) {
        if (apply_subsystem_cross_kick(bridge, &bridge->subsystems[i], dt_kick) != 0) {
            return -1;
        }
    }
    return 0;
}

static int bridge_integrate_all_to_target(struct reb_bridge* bridge, double target) {
    for (int i = 0; i < bridge->n_subsystems; i++) {
        enum REB_STATUS status = reb_simulation_integrate(bridge->subsystems[i].sub_sim, target);
        if (status != REB_STATUS_SUCCESS) return bridge_error("subsystem integration failed");
    }

    enum REB_STATUS main_status = reb_simulation_integrate(bridge->main_sim, target);
    if (main_status != REB_STATUS_SUCCESS) return bridge_error("main integration failed");
    return 0;
}

int reb_bridge_step(struct reb_bridge* bridge, double dt) {
    if (dt <= 0.0) return bridge_error("dt must be positive");

    if (reb_bridge_apply_cross_kick(bridge, 0.5 * dt) != 0) return -1;
    if (bridge_integrate_all_to_target(bridge, bridge->main_sim->t + dt) != 0) return -1;
    if (reb_bridge_apply_cross_kick(bridge, 0.5 * dt) != 0) return -1;
    return 0;
}

int reb_bridge_advance(struct reb_bridge* bridge, double duration) {
    if (duration < 0.0) return bridge_error("duration must be non-negative");
    if (duration == 0.0) return 0;

    const double target = bridge->main_sim->t + duration;
    double dt = target - bridge->main_sim->t;
    if (dt > bridge->dt_outer) dt = bridge->dt_outer;

    if (reb_bridge_apply_cross_kick(bridge, 0.5 * dt) != 0) return -1;

    /* Adjacent half-kicks are merged, preserving the KDK composition over many outer steps. */
    while (1) {
        const double step_target = bridge->main_sim->t + dt;
        if (bridge_integrate_all_to_target(bridge, step_target) != 0) return -1;

        if (bridge->main_sim->t >= target - 1.0e-14) break;

        double next_dt = target - bridge->main_sim->t;
        if (next_dt > bridge->dt_outer) next_dt = bridge->dt_outer;

        if (reb_bridge_apply_cross_kick(bridge, 0.5 * (dt + next_dt)) != 0) return -1;
        dt = next_dt;
    }

    if (reb_bridge_apply_cross_kick(bridge, 0.5 * dt) != 0) return -1;
    return 0;
}

int reb_bridge_snapshot(
    const struct reb_bridge* bridge,
    int subsystem_index,
    struct reb_bridge_snapshot* out
) {
    if (subsystem_index < 0 || subsystem_index >= bridge->n_subsystems) {
        return bridge_error("subsystem_index is out of bounds");
    }

    const struct reb_bridge_subsystem* sub = &bridge->subsystems[subsystem_index];
    const struct reb_simulation* main_sim = bridge->main_sim;
    const struct reb_simulation* sub_sim = sub->sub_sim;
    if (!valid_particle_index(main_sim, sub->host_index)) {
        return bridge_error("host_index is out of bounds during snapshot");
    }

    memset(out, 0, sizeof(*out));
    out->time = main_sim->t;
    out->host_index = sub->host_index;
    out->member_count = sub_sim->N;
    out->host = body_state(&main_sim->particles[sub->host_index]);

    if (main_sim->N > 0 && sub->host_index != 0) {
        out->host_primary_distance = particle_distance(
            &main_sim->particles[sub->host_index],
            &main_sim->particles[0]
        );
    }

    struct reb_vec3d bc;
    if (subsystem_barycenter(sub_sim, &bc) != 0) return -1;
    out->subsystem_barycenter_distance = sqrt(bc.x * bc.x + bc.y * bc.y + bc.z * bc.z);

    if (sub_sim->N > 0) {
        out->member0 = body_state(&sub_sim->particles[0]);
    }
    if (sub_sim->N > 1) {
        out->member1 = body_state(&sub_sim->particles[1]);
        out->subsystem_internal_distance = particle_distance(
            &sub_sim->particles[0],
            &sub_sim->particles[1]
        );
    }

    return 0;
}

struct reb_simulation* reb_bridge_make_main_sun_emb(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    reb_simulation_add_fmt(sim, "m", 1.0);
    reb_simulation_add_fmt(sim, "m a e", 3.0e-6, 1.0, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

struct reb_simulation* reb_bridge_make_sub_earth_moon(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }

    const double m_earth = 3.0e-6 * 0.987;
    const double m_moon = 3.0e-6 * 0.013;
    reb_simulation_add_fmt(sim, "m", m_earth);
    reb_simulation_add_fmt(sim, "m a e", m_moon, 0.00257, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

struct reb_bridge* reb_bridge_create_earth_moon(double dt_outer, int sub_ratio) {
    if (sub_ratio <= 0) {
        bridge_error("sub_ratio must be positive");
        return NULL;
    }

    struct reb_simulation* main_sim = reb_bridge_make_main_sun_emb();
    struct reb_simulation* sub_sim = reb_bridge_make_sub_earth_moon();
    if (!main_sim || !sub_sim) {
        if (main_sim) reb_simulation_free(main_sim);
        if (sub_sim) reb_simulation_free(sub_sim);
        return NULL;
    }

    struct reb_bridge* bridge = reb_bridge_create(main_sim, dt_outer);
    if (!bridge) {
        reb_simulation_free(main_sim);
        reb_simulation_free(sub_sim);
        return NULL;
    }

    reb_bridge_set_owns_main(bridge, 1);
    if (reb_bridge_add_subsystem(bridge, 1, sub_sim, dt_outer / (double)sub_ratio, 1) != 0) {
        reb_simulation_free(sub_sim);
        reb_bridge_free(bridge);
        return NULL;
    }
    return bridge;
}
