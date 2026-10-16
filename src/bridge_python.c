/* SPDX-License-Identifier: GPL-3.0-or-later */
/* The Python wrapper supplies symbols from its already-loaded REBOUND library.
 * No REBOUND source or second allocator is linked into this shared library. */
//this is  first warpper
#include <Python.h>
#include <stddef.h>
#include "bridge_internal.h"

typedef enum REB_STATUS (*integrate_fn)(struct reb_simulation*, double);
typedef void (*simulation_fn)(struct reb_simulation*);
typedef void (*transform_fn)(struct reb_particle*, const struct reb_particle*, const struct reb_particle*, unsigned int, unsigned int);

static integrate_fn loaded_integrate;
static simulation_fn loaded_synchronize;
static simulation_fn loaded_free;
static transform_fn loaded_transform;

REB_BRIDGE_API int reb_bridge_python_bind(void* integrate, void* synchronize, void* release, void* transform) {
    if (!integrate || !synchronize || !release || !transform) return -1;
    loaded_integrate = (integrate_fn)integrate;
    loaded_synchronize = (simulation_fn)synchronize;
    loaded_free = (simulation_fn)release;
    loaded_transform = (transform_fn)transform;
    return 0;
}

/* Keep every Python batch on the original absolute grid, including after a
 * checkpoint restart. The public C duration API retains relative semantics. */
REB_BRIDGE_API int reb_bridge_python_advance_steps(
    struct reb_bridge* bridge, double epoch,
    unsigned long long completed, unsigned long long count
) {
    return reb_bridge_internal_advance_grid(bridge, epoch, completed, count);
}

/* Only for validated checkpoints. Tiny accumulated barycentre drift is part
 * of the saved numerical state and must not be changed by recentering. */
REB_BRIDGE_API int reb_bridge_python_add_restored_subsystem(
    struct reb_bridge* bridge, int host_index, struct reb_simulation* sim, double dt
) {
    return reb_bridge_internal_add_subsystem(bridge, host_index, sim, dt, 0, 0);
}

/* setup.py renames these at compile time to avoid ELF symbol interposition. */
enum REB_STATUS reb_simulation_integrate(struct reb_simulation* sim, double target) {
    if (!loaded_integrate) return REB_STATUS_GENERIC_ERROR;
    if (sim->integrator == REB_INTEGRATOR_IAS15 || sim->integrator == REB_INTEGRATOR_BS) {
        /* Python accepts only autonomous Newtonian components, with no user
         * callbacks, variational equations or external ODEs. A local stage
         * clock lets adaptive solvers shrink dt below the ULP of an absolute
         * epoch without hanging or losing physical integration time. Their
         * derivative caches and timestep proposals do not contain an epoch.
         * This wrapper is not part of the independent C library. */
        const double epoch = sim->t;
        sim->t = 0.0;
        const enum REB_STATUS status = loaded_integrate(sim, target - epoch);
        /* 失败时也要保留实际走过的时间 */
        sim->t = status == REB_STATUS_SUCCESS ? target : epoch + sim->t;
        return status;
    }
    return loaded_integrate(sim, target);
}
void reb_simulation_synchronize(struct reb_simulation* sim) {
    if (loaded_synchronize) loaded_synchronize(sim);
}
void reb_simulation_free(struct reb_simulation* sim) {
    if (loaded_free) loaded_free(sim);
}
void reb_particles_transform_jacobi_to_inertial_posvel(
    struct reb_particle* inertial, const struct reb_particle* jacobi,
    const struct reb_particle* masses, unsigned int n, unsigned int n_active
) {
    if (loaded_transform) loaded_transform(inertial, jacobi, masses, n, n_active);
}

REB_BRIDGE_API const char* reb_bridge_python_rebound_version(void) { return "4.6.0"; }
REB_BRIDGE_API size_t reb_bridge_python_sizeof_simulation(void) { return sizeof(struct reb_simulation); }
REB_BRIDGE_API size_t reb_bridge_python_sizeof_particle(void) { return sizeof(struct reb_particle); }
REB_BRIDGE_API size_t reb_bridge_python_offset(const char* name) {
#define FIELD(field) if (strcmp(name, #field) == 0) return offsetof(struct reb_simulation, field)
    FIELD(t); FIELD(G); FIELD(dt); FIELD(dt_last_done); FIELD(N); FIELD(N_var);
    FIELD(N_active); FIELD(particles); FIELD(exact_finish_time); FIELD(integrator);
    FIELD(collision); FIELD(boundary); FIELD(ri_whfast); FIELD(ri_ias15);
    FIELD(additional_forces); FIELD(pre_timestep_modifications);
    FIELD(post_timestep_modifications); FIELD(heartbeat); FIELD(extras);
#undef FIELD
    return (size_t)-1;
}

static struct reb_bridge_whfast_cache* python_cache(
    struct reb_bridge* bridge, int index, struct reb_simulation** sim
) {
    if (!bridge || index < -1 || index >= bridge->n_subsystems) return NULL;
    if (index == -1) {
        *sim = bridge->main_sim;
        return &bridge->main_whfast_cache;
    }
    *sim = bridge->subsystems[index].sub_sim;
    return &bridge->subsystems[index].whfast_cache;
}

REB_BRIDGE_API int reb_bridge_python_cache_valid(struct reb_bridge* bridge, int index) {
    struct reb_simulation* sim;
    struct reb_bridge_whfast_cache* cache = python_cache(bridge, index, &sim);
    if (!cache) return -1;
    if (!cache->valid || cache->count != (int)sim->N ||
        cache->jacobi != sim->ri_whfast.p_jh || !sim->ri_whfast.p_jh ||
        sim->ri_whfast.N_allocated != sim->N || !sim->ri_whfast.is_synchronized ||
        sim->ri_whfast.recalculate_coordinates_this_timestep) return 0;
    for (unsigned int i = 0; i < sim->N; i++) {
        const struct reb_particle* p = &sim->particles[i];
        const struct reb_bridge_cartesian_state* s = &cache->state[i];
        if (p->m != s->m || p->x != s->x || p->y != s->y || p->z != s->z ||
            p->vx != s->vx || p->vy != s->vy || p->vz != s->vz) return 0;
    }
    return 1;
}

/* 这里只恢复缓存记录，不重算 Jacobi 坐标 */
REB_BRIDGE_API int reb_bridge_python_restore_cache(struct reb_bridge* bridge, int index, int valid) {
    struct reb_simulation* sim;
    struct reb_bridge_whfast_cache* cache = python_cache(bridge, index, &sim);
    if (!cache) return -1;
    cache->valid = 0;
    if (!valid) return 0;
    if (sim->integrator != REB_INTEGRATOR_WHFAST ||
        sim->ri_whfast.coordinates != REB_WHFAST_COORDINATES_JACOBI ||
        !sim->ri_whfast.p_jh || !sim->ri_whfast.is_synchronized ||
        sim->ri_whfast.recalculate_coordinates_this_timestep ||
        sim->ri_whfast.N_allocated != sim->N) return -1;
    if (cache->capacity < (int)sim->N) {
        struct reb_bridge_cartesian_state* state = realloc(cache->state, sizeof(*state) * sim->N);
        if (!state) return -1;
        cache->state = state;
        cache->capacity = (int)sim->N;
    }
    for (unsigned int i = 0; i < sim->N; i++) {
        const struct reb_particle* p = &sim->particles[i];
        struct reb_bridge_cartesian_state* s = &cache->state[i];
        s->m = p->m; s->x = p->x; s->y = p->y; s->z = p->z;
        s->vx = p->vx; s->vy = p->vy; s->vz = p->vz;
    }
    cache->count = (int)sim->N;
    cache->jacobi = sim->ri_whfast.p_jh;
    cache->valid = 1;
    return 0;
}

/* Exposed as a normal extension so wheels carry the platform/Python ABI tag.
 * Python loads its exported C entry points with ctypes, just as REBOUND does. */
static struct PyModuleDef module = {PyModuleDef_HEAD_INIT, "_bridge", NULL, -1, NULL};
PyMODINIT_FUNC PyInit__bridge(void) { return PyModule_Create(&module); }
