# Pluto-Charon Binary Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a Pluto-Charon validation suite that proves Bridge can model a true hierarchical binary whose host particle is the binary center of mass.

**Architecture:** Add an independent `validation/pluto_charon_binary/` validation area with one comparison executable, one fair WHFast speed benchmark executable, one sweep script, and one plotting script. The comparison executable constructs both a Bridge model and an IAS15 physical-particle reference from the same frozen synthetic initial conditions, reconstructs Pluto and Charon global states from the Bridge COM host plus subsystem states, and reports COM, absolute-body, relative-binary, conservation, and timing diagnostics.

**Tech Stack:** C99, CMake 3.20, vendored REBOUND, `rebound_bridge`, PowerShell validation runners, Python/matplotlib plotting.

---

## File Structure

- Create: `validation/pluto_charon_binary/compare_pluto_charon.c`
- Create: `validation/pluto_charon_binary/bench_pluto_charon.c`
- Create: `validation/pluto_charon_binary/run_pluto_charon_sweep.ps1`
- Create: `validation/pluto_charon_binary/plot_pluto_charon.py`
- Create: `validation/pluto_charon_binary/.gitignore`
- Create: `tests/pluto_charon_initial_state_test.cmake`
- Create: `tests/pluto_charon_csv_schema_test.cmake`
- Create: `tests/pluto_charon_bench_cli_test.cmake`
- Modify: `CMakeLists.txt`
- Modify: `README.md`

The new C files may duplicate small helper functions from `validation/solar_system_ias15/compare_solar_system_ias15.c`. Do not extract shared validation libraries in this implementation. The priority is explicit, inspectable coordinate handling for Pluto-Charon.

All Pluto-Charon initial conditions are frozen in source. The first version uses a synthetic but physically scaled coplanar circular setup:

```c
#define PC_M_SUN 1.0
#define PC_M_PLUTO 6.552776758009927e-09
#define PC_M_CHARON 7.975981533540863e-10
#define PC_M_BINARY (PC_M_PLUTO + PC_M_CHARON)
#define PC_R_COM_AU 39.482
#define PC_V_COM_AU_YR 0.9999546351918088
#define PC_A_BINARY_AU 0.00013099116924797247
#define PC_X_PLUTO_COM_AU (-1.421401166032323e-05)
#define PC_X_CHARON_COM_AU 0.00011677715758764924
#define PC_VY_PLUTO_COM_AU_YR (-0.005107263360501845)
#define PC_VY_CHARON_COM_AU_YR 0.041959420925182254
```

Outer-planet extension uses frozen circular coplanar perturbers:

```c
{ "jupiter", 9.547919e-4, 5.2044, 0.3 },
{ "saturn", 2.8588567e-4, 9.5826, 1.8 },
{ "uranus", 4.366244e-5, 19.2184, 3.1 },
{ "neptune", 5.151389e-5, 30.1104, 4.7 },
```

The Bridge main simulation must use a Pluto-Charon COM host with mass `PC_M_BINARY`. It must not use Pluto as the host.

## Task 1: Add the Initial-State Regression Test

**Files:**
- Create: `tests/pluto_charon_initial_state_test.cmake`
- Modify: `CMakeLists.txt`
- Create after RED: `validation/pluto_charon_binary/compare_pluto_charon.c`

- [ ] **Step 1: Write the failing CTest script**

Create `tests/pluto_charon_initial_state_test.cmake`:

```cmake
if(NOT DEFINED PLUTO_CHARON_EXE)
    message(FATAL_ERROR "PLUTO_CHARON_EXE is required")
endif()

if(NOT DEFINED OUT_CSV)
    message(FATAL_ERROR "OUT_CSV is required")
endif()

execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case minimal
        --initial-state-csv "${OUT_CSV}"
    RESULT_VARIABLE compare_result
    OUTPUT_VARIABLE compare_stdout
    ERROR_VARIABLE compare_stderr
)

if(NOT compare_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon initial state check failed: ${compare_stdout}\n${compare_stderr}")
endif()

foreach(metric IN ITEMS max_abs_mass_diff max_dr_au max_dv_au_per_yr max_com_dr_au max_com_dv_au_per_yr)
    if(NOT compare_stdout MATCHES "${metric}=([0-9.eE+-]+)")
        message(FATAL_ERROR "Missing ${metric} in output: ${compare_stdout}")
    endif()
    set("${metric}" "${CMAKE_MATCH_1}")
endforeach()

if(max_abs_mass_diff GREATER 1e-15)
    message(FATAL_ERROR "Initial masses differ: max_abs_mass_diff=${max_abs_mass_diff}")
endif()
if(max_dr_au GREATER 1e-13)
    message(FATAL_ERROR "Initial positions differ: max_dr_au=${max_dr_au}")
endif()
if(max_dv_au_per_yr GREATER 1e-13)
    message(FATAL_ERROR "Initial velocities differ: max_dv_au_per_yr=${max_dv_au_per_yr}")
endif()
if(max_com_dr_au GREATER 1e-13)
    message(FATAL_ERROR "Initial COM positions differ: max_com_dr_au=${max_com_dr_au}")
endif()
if(max_com_dv_au_per_yr GREATER 1e-13)
    message(FATAL_ERROR "Initial COM velocities differ: max_com_dv_au_per_yr=${max_com_dv_au_per_yr}")
endif()

file(READ "${OUT_CSV}" csv_text)
string(REGEX MATCH "^[^\r\n]*" csv_header "${csv_text}")

set(required_columns
    body
    bridge_m
    ias15_m
    mass_diff
    bridge_x
    bridge_y
    bridge_z
    bridge_vx
    bridge_vy
    bridge_vz
    ias15_x
    ias15_y
    ias15_z
    ias15_vx
    ias15_vy
    ias15_vz
    dr_au
    dv_au_per_yr
)

foreach(column IN LISTS required_columns)
    if(NOT ",${csv_header}," MATCHES ",${column},")
        message(FATAL_ERROR "CSV header is missing required column '${column}'. Header: ${csv_header}")
    endif()
endforeach()

foreach(body IN ITEMS pluto charon pc_com)
    if(NOT csv_text MATCHES "\n${body},")
        message(FATAL_ERROR "Initial-state CSV is missing row '${body}': ${csv_text}")
    endif()
endforeach()
```

- [ ] **Step 2: Wire the future executable and test into CMake**

In `CMakeLists.txt`, inside `if(REBOUND_BRIDGE_BUILD_VALIDATION)`, add:

```cmake
    add_executable(
        pluto_charon_compare
        validation/pluto_charon_binary/compare_pluto_charon.c
    )
    target_link_libraries(pluto_charon_compare PRIVATE rebound_bridge::rebound_bridge)
```

Inside `if(REBOUND_BRIDGE_BUILD_TESTS)`, add:

```cmake
    if(TARGET pluto_charon_compare)
        add_test(
            NAME pluto_charon_initial_state_test
            COMMAND
                ${CMAKE_COMMAND}
                -DPLUTO_CHARON_EXE=$<TARGET_FILE:pluto_charon_compare>
                -DOUT_CSV=${CMAKE_CURRENT_BINARY_DIR}/pluto_charon_initial_state_test.csv
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/pluto_charon_initial_state_test.cmake
        )
    endif()
```

- [ ] **Step 3: Run configure to verify RED**

Run:

```powershell
cmake -S . -B build-pluto-charon -DREBOUND_BRIDGE_BUILD_EXAMPLES=OFF -DREBOUND_BRIDGE_BUILD_VALIDATION=ON -DREBOUND_BRIDGE_BUILD_TESTS=ON -DREBOUND_BRIDGE_INSTALL=OFF
```

Expected:

```text
CMake configure fails because validation/pluto_charon_binary/compare_pluto_charon.c does not exist.
```

- [ ] **Step 4: Implement the minimal initial-state executable**

Create `validation/pluto_charon_binary/compare_pluto_charon.c` with these units:

```c
/* SPDX-License-Identifier: GPL-3.0-or-later */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rebound_bridge.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

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

enum pc_main_index {
    PC_MAIN_SUN = 0,
    PC_MAIN_PC_COM_MINIMAL = 1,
    PC_MAIN_PC_COM_OUTER = 5
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
```

Use helper functions:

```c
static double sqr(double x) { return x * x; }

static double distance_state(const struct pc_body_state* a, const struct pc_body_state* b) {
    return sqrt(sqr(a->x - b->x) + sqr(a->y - b->y) + sqr(a->z - b->z));
}

static double velocity_distance_state(const struct pc_body_state* a, const struct pc_body_state* b) {
    return sqrt(sqr(a->vx - b->vx) + sqr(a->vy - b->vy) + sqr(a->vz - b->vz));
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
```

Implement constructors:

```c
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

static struct reb_simulation* make_ias15_reference(enum pc_case_kind case_kind) {
    if (case_kind == PC_CASE_OUTER) return NULL;
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_ias15(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    reb_simulation_add_fmt(sim, "m", PC_M_SUN);
    add_pluto_charon_physical_global(sim);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_bridge_main(enum pc_case_kind case_kind, int* host_index) {
    if (case_kind == PC_CASE_OUTER) return NULL;
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->exact_finish_time = 1;
    reb_simulation_add_fmt(sim, "m", PC_M_SUN);
    *host_index = PC_MAIN_PC_COM_MINIMAL;
    add_cartesian(sim, PC_M_BINARY, PC_R_COM_AU, 0.0, 0.0, 0.0, PC_V_COM_AU_YR, 0.0);
    reb_simulation_move_to_com(sim);
    return sim;
}

static struct reb_simulation* make_binary_subsystem(void) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
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
```

Implement reconstruction:

```c
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
```

Implement `--initial-state-csv` mode so it writes rows for `pluto`, `charon`, and `pc_com`, prints all five metrics used by the test, and exits without integrating.

- [ ] **Step 5: Verify GREEN**

Run:

```powershell
cmake -S . -B build-pluto-charon -DREBOUND_BRIDGE_BUILD_EXAMPLES=OFF -DREBOUND_BRIDGE_BUILD_VALIDATION=ON -DREBOUND_BRIDGE_BUILD_TESTS=ON -DREBOUND_BRIDGE_INSTALL=OFF
cmake --build build-pluto-charon --config Release --target pluto_charon_compare
ctest --test-dir build-pluto-charon -C Release -R pluto_charon_initial_state_test --output-on-failure
```

Expected:

```text
100% tests passed, 0 tests failed out of 1
```

- [ ] **Step 6: Commit**

Run:

```powershell
git add CMakeLists.txt tests/pluto_charon_initial_state_test.cmake validation/pluto_charon_binary/compare_pluto_charon.c
git commit -m "test: add Pluto-Charon initial state validation"
```

Expected:

```text
Commit contains only the initial-state test, CMake target, and compare_pluto_charon.c.
```

## Task 2: Add Runtime CSV Diagnostics for the Minimal Case

**Files:**
- Create: `tests/pluto_charon_csv_schema_test.cmake`
- Modify: `CMakeLists.txt`
- Modify: `validation/pluto_charon_binary/compare_pluto_charon.c`

- [ ] **Step 1: Write the failing CSV schema test**

Create `tests/pluto_charon_csv_schema_test.cmake`:

```cmake
if(NOT DEFINED PLUTO_CHARON_EXE)
    message(FATAL_ERROR "PLUTO_CHARON_EXE is required")
endif()

if(NOT DEFINED OUT_CSV)
    message(FATAL_ERROR "OUT_CSV is required")
endif()

execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case minimal
        --years 1
        --samples 4
        --dt-outer-days 4
        --binary-ratio 64
        --bridge-method eos-lf4-2
        --csv "${OUT_CSV}"
    RESULT_VARIABLE compare_result
    OUTPUT_VARIABLE compare_stdout
    ERROR_VARIABLE compare_stderr
)

if(NOT compare_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon comparison failed: ${compare_stdout}\n${compare_stderr}")
endif()

foreach(metric IN ITEMS bridge_method bridge_total_seconds ias15_total_seconds overall_speedup max_pc_com_abs_error_au max_pc_relative_error_au)
    if(NOT compare_stdout MATCHES "${metric}=([^\\r\\n]+)")
        message(FATAL_ERROR "Missing ${metric} in output: ${compare_stdout}")
    endif()
endforeach()

file(READ "${OUT_CSV}" csv_text)
string(REGEX MATCH "^[^\r\n]*" csv_header "${csv_text}")

set(required_columns
    case_name
    bridge_method
    dt_outer_days
    binary_ratio
    time_yr
    bridge_step_seconds
    ias15_step_seconds
    step_speedup
    bridge_cumulative_seconds
    ias15_cumulative_seconds
    cumulative_speedup
    bridge_energy_rel
    ias15_energy_rel
    energy_rel_diff
    bridge_l_rel
    ias15_l_rel
    l_rel_diff
    bridge_p_rel
    ias15_p_rel
    p_rel_diff
    pc_com_abs_error_au
    pc_com_velocity_error_au_per_yr
    pluto_abs_error_au
    charon_abs_error_au
    pc_relative_error_au
    bridge_pc_distance_au
    ias15_pc_distance_au
    pc_distance_diff_au
    bridge_pc_a_au
    ias15_pc_a_au
    pc_a_diff_au
    bridge_pc_e
    ias15_pc_e
    pc_e_diff
    bridge_com_sun_a_au
    ias15_com_sun_a_au
    com_sun_a_diff_au
    bridge_com_sun_e
    ias15_com_sun_e
    com_sun_e_diff
)

foreach(column IN LISTS required_columns)
    if(NOT ",${csv_header}," MATCHES ",${column},")
        message(FATAL_ERROR "CSV header is missing required column '${column}'. Header: ${csv_header}")
    endif()
endforeach()

string(REGEX MATCHALL "\n" newline_matches "${csv_text}")
list(LENGTH newline_matches newline_count)
if(NOT newline_count EQUAL 5)
    message(FATAL_ERROR "Expected header plus 4 sample rows, found ${newline_count} newline(s)")
endif()
```

- [ ] **Step 2: Add the test to CMake**

Inside `if(TARGET pluto_charon_compare)` in the tests block:

```cmake
        add_test(
            NAME pluto_charon_csv_schema_test
            COMMAND
                ${CMAKE_COMMAND}
                -DPLUTO_CHARON_EXE=$<TARGET_FILE:pluto_charon_compare>
                -DOUT_CSV=${CMAKE_CURRENT_BINARY_DIR}/pluto_charon_schema_test.csv
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/pluto_charon_csv_schema_test.cmake
        )
```

- [ ] **Step 3: Run the test to verify RED**

Run:

```powershell
cmake --build build-pluto-charon --config Release --target pluto_charon_compare
ctest --test-dir build-pluto-charon -C Release -R pluto_charon_csv_schema_test --output-on-failure
```

Expected:

```text
Test fails because --csv runtime comparison mode is not implemented.
```

- [ ] **Step 4: Implement runtime integration and diagnostics**

In `compare_pluto_charon.c`, add:

```c
static int parse_double_value(const char* text, double* out);
static int parse_int_value(const char* text, int* out);
static int parse_case_kind(const char* text, enum pc_case_kind* out);
static int parse_bridge_method(const char* text, enum reb_bridge_method* out);
static const char* case_name(enum pc_case_kind case_kind);
static const char* bridge_method_name(enum reb_bridge_method method);
static double monotonic_seconds(void);
static struct pc_bridge_state make_bridge_state(enum pc_case_kind case_kind, double dt_outer_days, int binary_ratio);
static void free_bridge_state(struct pc_bridge_state* state);
static void extract_ias15_pair(const struct reb_simulation* sim, enum pc_case_kind case_kind, struct pc_pair_state* out);
static struct reb_orbit binary_orbit(const struct pc_pair_state* pair);
static struct reb_orbit com_sun_orbit(const struct pc_body_state* sun, const struct pc_body_state* com);
```

Use these index rules for IAS15 physical bodies:

```c
static int ias15_pluto_index(enum pc_case_kind case_kind) {
    return case_kind == PC_CASE_OUTER ? 5 : 1;
}

static int ias15_charon_index(enum pc_case_kind case_kind) {
    return case_kind == PC_CASE_OUTER ? 6 : 2;
}
```

Compute conservation from physical reconstructed bodies, not from the Bridge proxy host:

```c
struct pc_conservation {
    double energy;
    double lx, ly, lz;
    double px, py, pz;
};

static struct pc_conservation compute_conservation(const struct pc_body_state* bodies, int n_bodies) {
    struct pc_conservation out;
    memset(&out, 0, sizeof(out));
    for (int i = 0; i < n_bodies; i++) {
        const struct pc_body_state* a = &bodies[i];
        const double v2 = a->vx * a->vx + a->vy * a->vy + a->vz * a->vz;
        out.energy += 0.5 * a->m * v2;
        out.lx += a->m * (a->y * a->vz - a->z * a->vy);
        out.ly += a->m * (a->z * a->vx - a->x * a->vz);
        out.lz += a->m * (a->x * a->vy - a->y * a->vx);
        out.px += a->m * a->vx;
        out.py += a->m * a->vy;
        out.pz += a->m * a->vz;
        for (int j = i + 1; j < n_bodies; j++) {
            const struct pc_body_state* b = &bodies[j];
            const double dx = b->x - a->x;
            const double dy = b->y - a->y;
            const double dz = b->z - a->z;
            const double r = sqrt(dx * dx + dy * dy + dz * dz);
            out.energy -= REB_BRIDGE_G_AU_YR_MSUN * a->m * b->m / r;
        }
    }
    return out;
}
```

Write a CSV header exactly matching the schema test. Generate exactly `samples` rows at:

```c
const double target_time = years * (double)(sample_index + 1) / (double)samples;
```

For each sample:

1. Advance Bridge by `target_time - bridge.main_sim->t`.
2. Advance IAS15 to `target_time`.
3. Reconstruct Bridge Pluto, Charon, and COM.
4. Extract IAS15 Pluto, Charon, and mass-weighted COM.
5. Compute absolute Pluto error, absolute Charon error, COM position error, COM velocity error.
6. Compute relative vector error using `(charon - pluto)` in Bridge and IAS15.
7. Compute Pluto-Charon relative orbit elements using `reb_orbit_from_particle`.
8. Compute COM orbit around Sun using the reconstructed COM and physical Sun.
9. Compute conservation relative to the initial conservation values.
10. Track maxima for stdout.

Print this summary:

```text
case_name=<minimal|outer>
bridge_method=<order2|order4|eos-lf4-2>
years=<value>
samples=<value>
dt_outer_days=<value>
binary_ratio=<value>
binary_dt_days=<value>
bridge_total_seconds=<value>
ias15_total_seconds=<value>
overall_speedup=<value>
max_pc_com_abs_error_au=<value>
max_pc_com_velocity_error_au_per_yr=<value>
max_pluto_abs_error_au=<value>
max_charon_abs_error_au=<value>
max_pc_relative_error_au=<value>
max_abs_pc_distance_diff_au=<value>
max_abs_pc_a_diff_au=<value>
max_abs_pc_e_diff=<value>
max_abs_bridge_energy_rel=<value>
max_abs_bridge_l_rel=<value>
max_abs_bridge_p_rel=<value>
```

- [ ] **Step 5: Verify GREEN**

Run:

```powershell
cmake --build build-pluto-charon --config Release --target pluto_charon_compare
ctest --test-dir build-pluto-charon -C Release -R "pluto_charon_(initial_state|csv_schema)_test" --output-on-failure
```

Expected:

```text
100% tests passed, 0 tests failed out of 2
```

- [ ] **Step 6: Commit**

Run:

```powershell
git add CMakeLists.txt tests/pluto_charon_csv_schema_test.cmake validation/pluto_charon_binary/compare_pluto_charon.c
git commit -m "feat: add Pluto-Charon runtime diagnostics"
```

Expected:

```text
Commit contains minimal-case runtime diagnostics and CSV schema test.
```

## Task 3: Add the Outer-Planets Case

**Files:**
- Modify: `tests/pluto_charon_csv_schema_test.cmake`
- Modify: `validation/pluto_charon_binary/compare_pluto_charon.c`

- [ ] **Step 1: Extend the schema test to run the outer case**

In `tests/pluto_charon_csv_schema_test.cmake`, add a second `execute_process` block after the minimal case:

```cmake
set(OUTER_CSV "${OUT_CSV}.outer.csv")
execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case outer
        --years 1
        --samples 4
        --dt-outer-days 4
        --binary-ratio 64
        --bridge-method eos-lf4-2
        --csv "${OUTER_CSV}"
    RESULT_VARIABLE outer_result
    OUTPUT_VARIABLE outer_stdout
    ERROR_VARIABLE outer_stderr
)

if(NOT outer_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon outer comparison failed: ${outer_stdout}\n${outer_stderr}")
endif()

if(NOT outer_stdout MATCHES "case_name=outer")
    message(FATAL_ERROR "Outer run did not report case_name=outer: ${outer_stdout}")
endif()

file(READ "${OUTER_CSV}" outer_csv_text)
if(NOT outer_csv_text MATCHES "\nouter,eos-lf4-2,")
    message(FATAL_ERROR "Outer CSV is missing outer/eos-lf4-2 sample rows: ${outer_csv_text}")
endif()
```

- [ ] **Step 2: Run the test to verify RED**

Run:

```powershell
cmake --build build-pluto-charon --config Release --target pluto_charon_compare
ctest --test-dir build-pluto-charon -C Release -R pluto_charon_csv_schema_test --output-on-failure
```

Expected:

```text
Test fails because --case outer is parsed but outer planets are not inserted into both simulations.
```

- [ ] **Step 3: Implement frozen outer planets in both Bridge and IAS15**

In `compare_pluto_charon.c`, add:

```c
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
```

Call `add_outer_planets(sim)` in both `make_ias15_reference(PC_CASE_OUTER)` and `make_bridge_main(PC_CASE_OUTER)` before adding Pluto-Charon.

When computing conservation bodies for the outer case, include Sun, Jupiter, Saturn, Uranus, Neptune, reconstructed Pluto, and reconstructed Charon for Bridge. Include the same physical particles from IAS15 for the reference.

- [ ] **Step 4: Verify GREEN**

Run:

```powershell
cmake --build build-pluto-charon --config Release --target pluto_charon_compare
ctest --test-dir build-pluto-charon -C Release -R pluto_charon_csv_schema_test --output-on-failure
```

Expected:

```text
100% tests passed, 0 tests failed out of 1
```

- [ ] **Step 5: Commit**

Run:

```powershell
git add tests/pluto_charon_csv_schema_test.cmake validation/pluto_charon_binary/compare_pluto_charon.c
git commit -m "feat: add Pluto-Charon outer planet validation case"
```

Expected:

```text
Commit contains only outer-case support and the extended schema test.
```

## Task 4: Add Fair Pure-WHFast Speed Baseline

**Files:**
- Create: `tests/pluto_charon_bench_cli_test.cmake`
- Create: `validation/pluto_charon_binary/bench_pluto_charon.c`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Write the failing benchmark CLI test**

Create `tests/pluto_charon_bench_cli_test.cmake`:

```cmake
if(NOT DEFINED BENCH_EXE)
    message(FATAL_ERROR "BENCH_EXE is required")
endif()

execute_process(
    COMMAND "${BENCH_EXE}"
        --case minimal
        --years 0.01
        --dt-outer-days 4
        --binary-ratio 64
        --bridge-method eos-lf4-2
        --repeats 1
        --warmups 0
    RESULT_VARIABLE bench_result
    OUTPUT_VARIABLE bench_stdout
    ERROR_VARIABLE bench_stderr
)

if(NOT bench_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon benchmark failed: ${bench_stdout}\n${bench_stderr}")
endif()

foreach(metric IN ITEMS case_name bridge_method whfast_dt_days bridge_mean_seconds whfast_mean_seconds mean_speedup_whfast_over_bridge)
    if(NOT bench_stdout MATCHES "${metric}=([^\\r\\n]+)")
        message(FATAL_ERROR "Missing ${metric} in benchmark output: ${bench_stdout}")
    endif()
endforeach()

if(NOT bench_stdout MATCHES "case_name=minimal")
    message(FATAL_ERROR "Expected minimal benchmark output: ${bench_stdout}")
endif()
if(NOT bench_stdout MATCHES "bridge_method=eos-lf4-2")
    message(FATAL_ERROR "Expected eos-lf4-2 benchmark output: ${bench_stdout}")
endif()
if(NOT bench_stdout MATCHES "whfast_dt_days=0\\.062500000")
    message(FATAL_ERROR "Expected fair WHFast timestep 4/64 days: ${bench_stdout}")
endif()
```

- [ ] **Step 2: Wire the future benchmark target and test into CMake**

Inside `if(REBOUND_BRIDGE_BUILD_VALIDATION)`:

```cmake
    add_executable(
        pluto_charon_whfast_bridge_bench
        validation/pluto_charon_binary/bench_pluto_charon.c
    )
    target_link_libraries(pluto_charon_whfast_bridge_bench PRIVATE rebound_bridge::rebound_bridge)
```

Inside `if(REBOUND_BRIDGE_BUILD_TESTS)`:

```cmake
    if(TARGET pluto_charon_whfast_bridge_bench)
        add_test(
            NAME pluto_charon_bench_cli_test
            COMMAND
                ${CMAKE_COMMAND}
                -DBENCH_EXE=$<TARGET_FILE:pluto_charon_whfast_bridge_bench>
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/pluto_charon_bench_cli_test.cmake
        )
    endif()
```

- [ ] **Step 3: Run configure/build to verify RED**

Run:

```powershell
cmake -S . -B build-pluto-charon -DREBOUND_BRIDGE_BUILD_EXAMPLES=OFF -DREBOUND_BRIDGE_BUILD_VALIDATION=ON -DREBOUND_BRIDGE_BUILD_TESTS=ON -DREBOUND_BRIDGE_INSTALL=OFF
cmake --build build-pluto-charon --config Release --target pluto_charon_whfast_bridge_bench
```

Expected:

```text
Build fails because validation/pluto_charon_binary/bench_pluto_charon.c does not exist.
```

- [ ] **Step 4: Implement the benchmark executable**

Create `validation/pluto_charon_binary/bench_pluto_charon.c`. Copy the same frozen constants, case parser, method parser, `add_cartesian`, `add_outer_planets`, `make_bridge_main`, and `make_binary_subsystem` helpers from `compare_pluto_charon.c`.

Implement:

```c
static struct reb_simulation* make_full_whfast(enum pc_case_kind case_kind, double dt_days) {
    struct reb_simulation* sim = reb_simulation_create();
    if (!sim) return NULL;
    sim->G = REB_BRIDGE_G_AU_YR_MSUN;
    if (reb_bridge_set_integrator_whfast(sim) != 0) {
        reb_simulation_free(sim);
        return NULL;
    }
    sim->dt = dt_days / 365.25;
    sim->exact_finish_time = 1;
    reb_simulation_add_fmt(sim, "m", PC_M_SUN);
    if (case_kind == PC_CASE_OUTER) add_outer_planets(sim);
    add_pluto_charon_physical_global(sim);
    reb_simulation_move_to_com(sim);
    return sim;
}
```

Use this fair timestep:

```c
const double binary_dt_days = dt_outer_days / (double)binary_ratio;
const double whfast_dt_days = binary_dt_days;
```

Run repeats in this order:

1. Warmup pure WHFast.
2. Warmup Bridge.
3. Repeat pure WHFast.
4. Repeat Bridge.

Print per-repeat lines:

```text
repeat=<n> whfast_seconds=<value> bridge_seconds=<value> speedup_whfast_over_bridge=<value>
```

Print summary lines:

```text
case_name=<minimal|outer>
years=<value>
dt_outer_days=<value>
binary_ratio=<value>
bridge_method=<order2|order4|eos-lf4-2>
whfast_dt_days=<value>
bridge_mean_seconds=<value>
whfast_mean_seconds=<value>
mean_speedup_whfast_over_bridge=<value>
bridge_best_seconds=<value>
whfast_best_seconds=<value>
best_speedup_whfast_over_bridge=<value>
```

- [ ] **Step 5: Verify GREEN**

Run:

```powershell
cmake --build build-pluto-charon --config Release --target pluto_charon_whfast_bridge_bench
ctest --test-dir build-pluto-charon -C Release -R pluto_charon_bench_cli_test --output-on-failure
```

Expected:

```text
100% tests passed, 0 tests failed out of 1
```

- [ ] **Step 6: Commit**

Run:

```powershell
git add CMakeLists.txt tests/pluto_charon_bench_cli_test.cmake validation/pluto_charon_binary/bench_pluto_charon.c
git commit -m "feat: add Pluto-Charon WHFast speed baseline"
```

Expected:

```text
Commit contains only the benchmark target, benchmark CLI test, and benchmark source.
```

## Task 5: Add Sweep Runner and Plotting

**Files:**
- Create: `validation/pluto_charon_binary/run_pluto_charon_sweep.ps1`
- Create: `validation/pluto_charon_binary/plot_pluto_charon.py`
- Create: `validation/pluto_charon_binary/.gitignore`

- [ ] **Step 1: Write the runner script**

Create `validation/pluto_charon_binary/run_pluto_charon_sweep.ps1`:

```powershell
param(
    [double]$Years = 100.0,
    [int]$Samples = 200,
    [int]$BinaryRatio = 64,
    [double[]]$DtOuterDays = @(2.0, 4.0, 8.0, 16.0),
    [string[]]$Cases = @("minimal", "outer"),
    [string[]]$Methods = @("order2", "eos-lf4-2"),
    [int]$CaseTimeoutSeconds = 600,
    [string]$BuildDir = "build-pluto-charon-validation",
    [string]$OutDir = "validation/pluto_charon_binary/benchmark_figures"
)

$ErrorActionPreference = "Stop"
$culture = [System.Globalization.CultureInfo]::InvariantCulture

function Format-Number([double]$value) {
    return $value.ToString("R", $culture)
}

function Read-KeyValueOutput([string]$text) {
    $values = @{}
    foreach ($line in ($text -split "`r?`n")) {
        if ($line -match "^([^=]+)=(.*)$") {
            $values[$matches[1]] = $matches[2]
        }
    }
    return $values
}

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$caseDir = Join-Path $OutDir "cases"
New-Item -ItemType Directory -Force -Path $caseDir | Out-Null

cmake -S . -B $BuildDir -DREBOUND_BRIDGE_BUILD_EXAMPLES=OFF -DREBOUND_BRIDGE_BUILD_VALIDATION=ON -DREBOUND_BRIDGE_BUILD_TESTS=OFF -DREBOUND_BRIDGE_INSTALL=OFF | Write-Host
cmake --build $BuildDir --config Release --target pluto_charon_compare pluto_charon_whfast_bridge_bench | Write-Host

$compareExe = Join-Path $BuildDir "Release/pluto_charon_compare.exe"
$benchExe = Join-Path $BuildDir "Release/pluto_charon_whfast_bridge_bench.exe"
$summaryPath = Join-Path $OutDir "pluto_charon_sweep.csv"
$benchPath = Join-Path $OutDir "pluto_charon_whfast_bridge_speed.csv"

$summaryHeader = @(
    "case",
    "status",
    "case_name",
    "bridge_method",
    "years",
    "samples",
    "dt_outer_days",
    "binary_ratio",
    "binary_dt_days",
    "bridge_seconds",
    "ias15_seconds",
    "ias15_over_bridge_speedup",
    "max_pc_com_abs_error_au",
    "max_pc_com_velocity_error_au_per_yr",
    "max_pluto_abs_error_au",
    "max_charon_abs_error_au",
    "max_pc_relative_error_au",
    "max_abs_pc_distance_diff_au",
    "max_abs_pc_a_diff_au",
    "max_abs_pc_e_diff",
    "max_abs_bridge_energy_rel",
    "max_abs_bridge_l_rel",
    "max_abs_bridge_p_rel",
    "case_csv"
)

Set-Content -LiteralPath $summaryPath -Value (($summaryHeader | ForEach-Object { '"' + $_ + '"' }) -join ",")
Set-Content -LiteralPath $benchPath -Value '"case","case_name","bridge_method","years","dt_outer_days","binary_ratio","whfast_dt_days","bridge_mean_seconds","whfast_mean_seconds","mean_speedup_whfast_over_bridge","bridge_best_seconds","whfast_best_seconds","best_speedup_whfast_over_bridge"'

foreach ($caseName in $Cases) {
    foreach ($method in $Methods) {
        foreach ($dt in $DtOuterDays) {
            $caseId = "{0}_{1}_dt{2}_r{3}" -f $caseName, ($method -replace "-", "_"), (($dt.ToString("0.###", $culture)) -replace "\.", "p"), $BinaryRatio
            $csvPath = Join-Path $caseDir "$caseId.csv"
            $args = @(
                "--case", $caseName,
                "--years", (Format-Number $Years),
                "--samples", $Samples.ToString($culture),
                "--dt-outer-days", (Format-Number $dt),
                "--binary-ratio", $BinaryRatio.ToString($culture),
                "--bridge-method", $method,
                "--csv", $csvPath
            )

            $stdout = ""
            $stderr = ""
            $status = "complete"
            try {
                $job = Start-Job -ScriptBlock {
                    param($exe, $argv)
                    & $exe @argv
                } -ArgumentList $compareExe, $args
                if (Wait-Job $job -Timeout $CaseTimeoutSeconds) {
                    $stdout = Receive-Job $job | Out-String
                    if ($job.State -ne "Completed") { $status = "failed" }
                } else {
                    Stop-Job $job
                    $status = "timeout"
                    $stdout = ""
                }
                Remove-Job $job -Force
            } catch {
                $status = "failed"
                $stderr = $_.Exception.Message
            }

            $values = Read-KeyValueOutput $stdout
            $row = @(
                $caseId,
                $status,
                $values["case_name"],
                $values["bridge_method"],
                $values["years"],
                $values["samples"],
                $values["dt_outer_days"],
                $values["binary_ratio"],
                $values["binary_dt_days"],
                $values["bridge_total_seconds"],
                $values["ias15_total_seconds"],
                $values["overall_speedup"],
                $values["max_pc_com_abs_error_au"],
                $values["max_pc_com_velocity_error_au_per_yr"],
                $values["max_pluto_abs_error_au"],
                $values["max_charon_abs_error_au"],
                $values["max_pc_relative_error_au"],
                $values["max_abs_pc_distance_diff_au"],
                $values["max_abs_pc_a_diff_au"],
                $values["max_abs_pc_e_diff"],
                $values["max_abs_bridge_energy_rel"],
                $values["max_abs_bridge_l_rel"],
                $values["max_abs_bridge_p_rel"],
                $csvPath
            )
            Add-Content -LiteralPath $summaryPath -Value (($row | ForEach-Object { '"' + ($_ -replace '"', '""') + '"' }) -join ",")

            $benchArgs = @(
                "--case", $caseName,
                "--years", (Format-Number $Years),
                "--dt-outer-days", (Format-Number $dt),
                "--binary-ratio", $BinaryRatio.ToString($culture),
                "--bridge-method", $method,
                "--repeats", "3",
                "--warmups", "1"
            )
            $benchOut = & $benchExe @benchArgs | Out-String
            $bench = Read-KeyValueOutput $benchOut
            $benchRow = @(
                $caseId,
                $bench["case_name"],
                $bench["bridge_method"],
                $bench["years"],
                $bench["dt_outer_days"],
                $bench["binary_ratio"],
                $bench["whfast_dt_days"],
                $bench["bridge_mean_seconds"],
                $bench["whfast_mean_seconds"],
                $bench["mean_speedup_whfast_over_bridge"],
                $bench["bridge_best_seconds"],
                $bench["whfast_best_seconds"],
                $bench["best_speedup_whfast_over_bridge"]
            )
            Add-Content -LiteralPath $benchPath -Value (($benchRow | ForEach-Object { '"' + ($_ -replace '"', '""') + '"' }) -join ",")
        }
    }
}

python validation/pluto_charon_binary/plot_pluto_charon.py --summary-csv $summaryPath --speed-csv $benchPath --outdir $OutDir
```

- [ ] **Step 2: Write the plotting script**

Create `validation/pluto_charon_binary/plot_pluto_charon.py` with command-line options:

```python
parser.add_argument("--summary-csv", required=True)
parser.add_argument("--speed-csv", required=True)
parser.add_argument("--outdir", required=True)
```

Generate these files:

```text
pluto_charon_com_vs_relative_error.png
pluto_charon_orbital_elements.png
pluto_charon_conservation.png
pluto_charon_speed_tradeoff.png
pluto_charon_dt_convergence.png
```

Each plot must separate `minimal` and `outer` cases and distinguish `order2` from `eos-lf4-2`.

Use `csv.DictReader`, `matplotlib.pyplot`, and semilog-y axes for absolute error plots. Use only columns listed in Task 2 and Task 5.

- [ ] **Step 3: Add output ignores**

Create `validation/pluto_charon_binary/.gitignore`:

```gitignore
/out/
/benchmark_figures/cases/
*.pyc
__pycache__/
```

- [ ] **Step 4: Run a short sweep**

Run:

```powershell
.\validation\pluto_charon_binary\run_pluto_charon_sweep.ps1 -Years 1 -Samples 8 -DtOuterDays @(4.0) -Cases @("minimal") -Methods @("order2","eos-lf4-2") -CaseTimeoutSeconds 120
```

Expected:

```text
validation/pluto_charon_binary/benchmark_figures/pluto_charon_sweep.csv exists
validation/pluto_charon_binary/benchmark_figures/pluto_charon_whfast_bridge_speed.csv exists
validation/pluto_charon_binary/benchmark_figures/pluto_charon_com_vs_relative_error.png exists
validation/pluto_charon_binary/benchmark_figures/pluto_charon_speed_tradeoff.png exists
```

- [ ] **Step 5: Commit**

Run:

```powershell
git add validation/pluto_charon_binary/run_pluto_charon_sweep.ps1 validation/pluto_charon_binary/plot_pluto_charon.py validation/pluto_charon_binary/.gitignore
git commit -m "feat: add Pluto-Charon validation sweep"
```

Expected:

```text
Commit contains only runner, plotter, and validation output ignore rules.
```

## Task 6: Document and Verify the Validation Suite

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Document the Pluto-Charon validation**

Add a `Pluto-Charon binary validation` subsection to `README.md`:

````markdown
### Pluto-Charon Binary Validation

`validation/pluto_charon_binary` validates Bridge on a true binary subsystem.
The main simulation uses a Pluto-Charon center-of-mass host particle, while the
subsystem resolves Pluto and Charon in their own COM frame.

The key diagnostic is the separation between:

- Pluto-Charon COM error, which diagnoses external coupling and cross-kick behavior.
- Pluto-Charon relative-orbit error, which diagnoses internal binary drift and splitting error.

Run a short validation sweep:

```powershell
.\validation\pluto_charon_binary\run_pluto_charon_sweep.ps1 -Years 1 -Samples 8 -DtOuterDays @(4.0) -Cases @("minimal") -Methods @("order2","eos-lf4-2")
```

Run the CTest coverage:

```powershell
cmake -S . -B build-pluto-charon -DREBOUND_BRIDGE_BUILD_EXAMPLES=OFF -DREBOUND_BRIDGE_BUILD_VALIDATION=ON -DREBOUND_BRIDGE_BUILD_TESTS=ON -DREBOUND_BRIDGE_INSTALL=OFF
cmake --build build-pluto-charon --config Release
ctest --test-dir build-pluto-charon -C Release -R pluto_charon --output-on-failure
```
````

- [ ] **Step 2: Run full targeted verification**

Run:

```powershell
cmake -S . -B build-pluto-charon -DREBOUND_BRIDGE_BUILD_EXAMPLES=OFF -DREBOUND_BRIDGE_BUILD_VALIDATION=ON -DREBOUND_BRIDGE_BUILD_TESTS=ON -DREBOUND_BRIDGE_INSTALL=OFF
cmake --build build-pluto-charon --config Release
ctest --test-dir build-pluto-charon -C Release -R pluto_charon --output-on-failure
git diff --check
```

Expected:

```text
All pluto_charon tests pass.
git diff --check reports no whitespace errors. CRLF warnings are acceptable on this Windows checkout.
```

- [ ] **Step 3: Run the short smoke sweep**

Run:

```powershell
.\validation\pluto_charon_binary\run_pluto_charon_sweep.ps1 -Years 1 -Samples 8 -DtOuterDays @(4.0) -Cases @("minimal") -Methods @("order2","eos-lf4-2") -CaseTimeoutSeconds 120
```

Expected:

```text
Both minimal cases complete.
Summary CSV has two complete rows.
At least five PNG files exist in validation/pluto_charon_binary/benchmark_figures.
```

- [ ] **Step 4: Remove temporary build directory**

Run:

```powershell
$target = Join-Path (Get-Location) 'build-pluto-charon'
$resolved = Resolve-Path -LiteralPath $target -ErrorAction Stop
if ($resolved.Path -ne $target) { throw "Unexpected build path: $($resolved.Path)" }
Remove-Item -LiteralPath $resolved.Path -Recurse -Force
```

Expected:

```text
build-pluto-charon no longer exists.
```

- [ ] **Step 5: Commit documentation**

Run:

```powershell
git add README.md
git commit -m "docs: describe Pluto-Charon binary validation"
```

Expected:

```text
Commit contains only README.md.
```

## Final Verification Checklist

- [ ] `ctest --test-dir build-pluto-charon -C Release -R pluto_charon --output-on-failure` passes before build directory cleanup.
- [ ] `validation/pluto_charon_binary/benchmark_figures/pluto_charon_sweep.csv` contains minimal `order2` and `eos-lf4-2` rows from the smoke sweep.
- [ ] `validation/pluto_charon_binary/benchmark_figures/pluto_charon_whfast_bridge_speed.csv` contains fair pure-WHFast speed rows.
- [ ] `validation/pluto_charon_binary/benchmark_figures/pluto_charon_com_vs_relative_error.png` exists.
- [ ] `git diff --check` has no whitespace errors.
- [ ] No `pluto_charon_compare` or `pluto_charon_whfast_bridge_bench` process remains running.
- [ ] Temporary build directories created by this plan are removed.
