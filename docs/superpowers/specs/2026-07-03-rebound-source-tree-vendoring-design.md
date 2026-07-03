# REBOUND Source-tree Vendoring Design

Date: 2026-07-03
Branch: whfasttest

## Goal

Remove the build-time dependency on an external REBOUND checkout or CMake
FetchContent download by vendoring the complete REBOUND upstream `src/` source
tree into this repository.

This is not a WHFast extraction. The implementation will embed REBOUND's source
tree as a local third-party source tree and compile it from this repository.
`rebound_bridge` will keep using the official REBOUND C API, including
`struct reb_simulation`, `struct reb_particle`, `reb_simulation_add_fmt()`,
WHFast, and IAS15.

## Core Decision

Use source-tree vendoring:

- Copy the entire upstream REBOUND `src/` directory into this repository.
- Preserve upstream source files exactly unless a build-system integration
  change is unavoidable.
- Preserve all upstream comments verbatim. Do not rewrite, translate, reformat,
  summarize, or delete copied upstream comments.
- Keep upstream filenames and internal organization intact.
- Stop trying to split WHFast into `whfast.c`, `whfast.h`, and
  `Keplersolve.c`. REBOUND's WHFast file stays as upstream provides it.

The vendored tree becomes the fixed source used by this project. `_deps/` can
remain a local cache but is no longer part of normal configure/build behavior.

## Source Organization

The project will add a vendored REBOUND directory:

```text
third_party/
  rebound/
    src/
      rebound.h
      rebound.c
      integrator_whfast.c
      integrator_whfast.h
      integrator_ias15.c
      integrator_ias15.h
      particle.c
      particle.h
      tools.c
      tools.h
      ...
```

The directory is intentionally a full source-tree copy. It is not a curated
WHFast subset.

The bridge source remains:

```text
include/
  rebound_bridge.h
src/
  rebound_bridge.c
```

`include/rebound_bridge.h` will continue including `rebound.h`, and downstream
code will continue interacting with REBOUND simulation objects directly.

## Build Changes

CMake will stop accepting or fetching an external `REBOUND_SRC_DIR` for the
default build.

Removed default behavior:

- No `FetchContent` download of REBOUND.
- No required `REBOUND_SRC_DIR`.
- No dependency on `_deps/rebound-upstream`.
- No configure-time validation against an external `rebound.h`.

New default behavior:

- Set the REBOUND source directory internally to
  `third_party/rebound/src`.
- Build `rebound_bridge_rebound` from the vendored REBOUND source files.
- Build `rebound_bridge` and link it against `rebound_bridge_rebound`.
- Install the vendored REBOUND headers from `third_party/rebound/src`.

The first implementation may keep the existing target names:

- `rebound_bridge::rebound`
- `rebound_bridge::rebound_bridge`

This minimizes disruption while removing the external source-tree dependency.

## API Impact

There is no intended public API break in this design.

The following existing APIs remain valid:

- `reb_bridge_create(struct reb_simulation* main_sim, double dt_outer)`
- `reb_bridge_add_subsystem(..., struct reb_simulation* sub_sim, ...)`
- `reb_bridge_main_sim()`
- `reb_bridge_sub_sim()`
- `reb_bridge_set_integrator_whfast()`
- `reb_bridge_set_integrator_ias15()`
- REBOUND APIs such as `reb_simulation_create()`,
  `reb_simulation_add_fmt()`, `reb_simulation_integrate()`,
  `reb_simulation_move_to_com()`, and `reb_orbit_from_particle()`

Examples and validation programs keep the same REBOUND-facing API because they
continue to use official REBOUND types. Any required changes are limited to
build paths or include paths.

## Vendoring Rules

Copied REBOUND files are third-party source.

Rules for copied upstream files:

- Preserve upstream comments verbatim.
- Preserve copyright and license headers verbatim.
- Do not change comment wording, formatting, language, or ordering.
- Do not split upstream implementation files as part of this vendoring pass.
- Do not rename upstream functions or structs.
- Do not prune upstream files in the first pass.

If build integration later requires a local patch to an upstream file, the patch
must be minimal and documented separately from normal bridge code changes.

## Validation

Baseline validation after implementation:

- Configure a clean build without setting `REBOUND_SRC_DIR`.
- Confirm CMake does not download REBOUND.
- Confirm CMake does not read `_deps/rebound-upstream`.
- Build `rebound_bridge`, examples, and validation targets that previously
  worked.
- Run the Earth-Moon example.
- Run the Jovian Laplace example.
- Run the solar-system bridge example.
- Run the existing WHFast and IAS15 validation tools if they build with the
  vendored tree.

Success condition: the project builds and runs from repository-contained source
only, while preserving the existing REBOUND-facing bridge API.

## Migration Notes

This design supersedes the earlier standalone-WHFast extraction idea.

The previous plan to create `src/integrator/WHFast/whfast.c`,
`src/integrator/WHFast/Keplersolve.c`, and a local
`reb_bridge_simulation` API is intentionally abandoned. Whole source-tree
vendoring is the simpler and more reliable dependency-removal strategy for the
current codebase.
