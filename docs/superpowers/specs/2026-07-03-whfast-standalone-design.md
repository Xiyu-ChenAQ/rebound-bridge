# Standalone WHFast Design

Date: 2026-07-03
Branch: whfasttest

## Goal

Remove the build-time dependency on the REBOUND 5 upstream source tree while
preserving the REBOUND-style programming model needed by `rebound_bridge`.

The implementation will vendor the relevant REBOUND source code into this
repository, preserving upstream comments and copyright/license headers in
copied code. The bridge will compile only the local files it needs instead of
globbing and building the entire REBOUND source directory.

## Scope

In scope:

- Add a local simulation/particle layer compatible with the bridge's current
  needs.
- Add a local integrator layout:
  - `src/integrator/WHFast/whfast.c`
  - `src/integrator/WHFast/whfast.h`
  - `src/integrator/WHFast/Keplersolve.c`
  - `src/integrator/WHFast/Keplersolve.h`
  - `src/integrator/IAS15/ias15.c`
  - `src/integrator/IAS15/ias15.h`
- Copy REBOUND WHFast implementation code directly and retain all upstream
  comments in copied blocks.
- Split WHFast's Kepler solver into `Keplersolve.c`.
- Remove CMake `FetchContent`, `REBOUND_SRC_DIR`, and the separate
  `rebound_bridge_rebound` target.
- Keep examples source-compatible where practical by providing
  `reb_bridge_simulation_add_fmt(..., "m a e", ...)` and Cartesian particle
  insertion helpers.

Out of scope for the first implementation pass:

- Full REBOUND feature compatibility.
- MPI, collisions, boundaries, visualization, archives, heartbeat callbacks,
  non-Newtonian forces, and tree gravity.
- WHFast512.
- A working standalone IAS15 integrator. The IAS15 folder exists as the future
  location, but validation should not depend on it in this pass.

## Source Organization

The repository will move away from depending on `_deps/rebound-upstream/src`.
The local runtime will be organized as:

```text
include/
  rebound_bridge.h

src/
  rebound_bridge.c
  simulation.c
  simulation.h
  particle.c
  particle.h
  tools.c
  tools.h
  integrator/
    WHFast/
      whfast.c
      whfast.h
      Keplersolve.c
      Keplersolve.h
    IAS15/
      ias15.c
      ias15.h
```

`WHFast` and `IAS15` intentionally use directory names matching the official
integrator names. File contents may preserve REBOUND-style function naming
internally, but public bridge-facing APIs will remain under the
`reb_bridge_` prefix to avoid implying ABI compatibility with upstream
REBOUND.

## Public API Shape

`include/rebound_bridge.h` will stop including `rebound.h`.

The public types become local:

```c
struct reb_bridge_particle {
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    double m;
};

struct reb_bridge_simulation;
```

Bridge APIs that currently accept `struct reb_simulation*` will accept
`struct reb_bridge_simulation*` instead:

- `reb_bridge_create()`
- `reb_bridge_add_subsystem()`
- `reb_bridge_main_sim()`
- `reb_bridge_sub_sim()`

Local simulation helpers will provide the REBOUND-like operations currently
used by examples and bridge internals:

- `reb_bridge_simulation_create()`
- `reb_bridge_simulation_free()`
- `reb_bridge_simulation_add_fmt()`
- `reb_bridge_simulation_integrate()`
- `reb_bridge_simulation_synchronize()`
- `reb_bridge_simulation_move_to_com()`

The old helper functions `reb_bridge_set_integrator_whfast()` and
`reb_bridge_set_integrator_ias15()` will be removed. WHFast is the built-in
first-pass integrator, so callers do not need to select it.

## `add_fmt` Compatibility

The local `reb_bridge_simulation_add_fmt()` will follow the parameter-selection
style of REBOUND's `reb_simulation_add_fmt()`.

First-pass supported formats:

- `"m"`
- `"m x y z vx vy vz"`
- `"m a e"`

For `"m a e"`, missing orbital orientation parameters default to zero:

- `inc = 0`
- `Omega = 0`
- `omega = 0`
- `f = 0`

The conversion uses the simulation's first particle as the primary when no
explicit primary is provided, matching the usage pattern in the current
examples. Later extensions can add `inc`, `Omega`, `omega`, `f`, `M`, `E`,
`primary`, and other upstream tokens.

Unsupported tokens should fail clearly instead of being silently ignored.

## WHFast Porting Rules

The WHFast implementation will be copied from REBOUND source code as directly
as possible:

- Preserve upstream comments in copied code blocks.
- Preserve upstream copyright/license headers in copied files or copied
  sections.
- Keep algorithmic structure recognizable to ease future comparison with
  REBOUND.
- Adapt type names and includes only where needed to target
  `reb_bridge_simulation` and `reb_bridge_particle`.
- Do not rewrite the copied WHFast algorithm into a different implementation
  during the dependency-removal pass.

The current upstream `integrator_whfast.c` contains the Kepler solver inline.
This project will split that logic into:

- `Keplersolve.c`: universal-variable Kepler solver and directly required
  helper functions/tables.
- `whfast.c`: WHFast initialization, coordinate transforms, interaction steps,
  jump steps, Kepler-step dispatch, synchronization, and integration loop.

The first implementation will prioritize the subset needed by the bridge's
current fixed-timestep Newtonian integrations. Unused upstream features can be
left out only if doing so does not require rewriting copied WHFast logic into a
different algorithm. If a feature is omitted, the omission must be explicit in
code comments or documentation.

## Build Changes

CMake will build a single local library target:

- `rebound_bridge`
- alias: `rebound_bridge::rebound_bridge`

Removed build behavior:

- No `FetchContent` for REBOUND.
- No `REBOUND_SRC_DIR`.
- No glob over upstream REBOUND `.c` files.
- No `rebound_bridge::rebound` target.
- No install step that copies upstream REBOUND headers.

The build will explicitly list local source files so accidental extra source
files do not change the library contents.

## Validation

Baseline validation after implementation:

- Configure and build with examples enabled.
- Run the Earth-Moon example.
- Run the solar-system bridge example if it still compiles under the supported
  `add_fmt` subset.
- Confirm the library no longer compiles or links upstream REBOUND sources.

IAS15 comparison tools can remain in the tree but will not be part of the
default validation build until local IAS15 exists or a legacy REBOUND-optional
validation mode is added.

## Migration Notes

Existing code using raw REBOUND types will need source changes:

- `struct reb_simulation` becomes `struct reb_bridge_simulation`.
- `struct reb_particle` becomes `struct reb_bridge_particle`.
- `reb_simulation_*` calls become `reb_bridge_simulation_*` calls.

The examples will be updated during implementation to use the local API.
This is an intentional API break to remove the upstream dependency cleanly.
