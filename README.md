# rebound_bridge

An operator-splitting scheduler for coupling independent REBOUND simulations.

It only depends on the REBOUND C API and can be moved into its own repository or built as a
separate static/shared library later.

## Model

`rebound_bridge` coordinates one main REBOUND simulation and one or more
subsystem simulations.

- The main simulation contains slow bodies and one host particle per subsystem.
- Each subsystem simulation contains the resolved internal bodies around its
  own barycenter.
- The bridge applies KDK cross kicks, then advances every REBOUND simulation
  to the same global time.

The Earth-Moon helpers are examples only.  The core API works with any
`host_index + sub_sim` pair.

## Using REBOUND

`rebound_bridge` is built against the REBOUND 5 C source tree.  By default,
CMake downloads REBOUND 5.0.0 from the upstream REBOUND repository and compiles
it as part of the build.

The default build is:

```powershell
cmake -S . -B build
cmake --build build --config Release
```

If you already have a local REBOUND 5 source tree, pass the directory that
contains `rebound.h` and the REBOUND `.c` files:

```powershell
cmake -S . -B build -DREBOUND_SRC_DIR="C:\path\to\rebound\src"
cmake --build build --config Release
```

If this repository contains `_deps/rebound-upstream/src/rebound.h` and that
copy provides the REBOUND 5 API, CMake uses it automatically.  Older REBOUND 4
source trees are rejected.

The build creates two static targets:

- `rebound_bridge::rebound` compiles the REBOUND C sources from
  `REBOUND_SRC_DIR`.
- `rebound_bridge::rebound_bridge` compiles this bridge library and links it
  against `rebound_bridge::rebound`.

Downstream CMake projects should link against:

```cmake
target_link_libraries(my_target PRIVATE rebound_bridge::rebound_bridge)
```

Manual builds need both include paths and both source sets:

```powershell
cl /Iinclude /IC:\path\to\rebound\src ^
  examples\solar_system_bridge.c src\rebound_bridge.c C:\path\to\rebound\src\*.c
```

At runtime, users still create normal REBOUND `struct reb_simulation*`
instances and choose REBOUND integrators as usual.  The bridge stores pointers
to those simulations and calls `reb_simulation_integrate()` during each bridge
step.

`dt_outer` and `dt_inner` are also written into REBOUND simulations:

- `reb_bridge_create(main_sim, dt_outer)` sets `main_sim->dt = dt_outer`.
- `reb_bridge_add_subsystem(..., sub_sim, dt_inner, ...)` sets
  `sub_sim->dt = dt_inner`.

For fixed-timestep REBOUND integrators such as WHFast, these values are the
actual REBOUND timesteps used during the drift integrations.  The bridge then
uses `dt_outer` as the interval between cross-system kick exchanges.  Therefore
`dt_outer` is both the main simulation timestep and the bridge coupling
timestep, while `dt_inner` is the subsystem timestep.  The current API requires
`dt_outer / dt_inner` to be an integer.

## Validation Notes

`rebound_bridge` is a scheduler around REBOUND simulations, not a standalone
replacement integrator.  A typical configuration uses WHFast inside the main
simulation and inside each subsystem, while the bridge supplies the symmetric
KDK coupling between subsystem barycenters and resolved subsystem bodies.

The validation target is long-term structural fidelity rather than pointwise
agreement with a high-order non-symplectic reference at every output time.
In practice this means checking that resolved subsystems stay bound, resonant
structure is preserved, and conserved quantities remain bounded over long
integrations.

Current validated setups:

- Forward-then-backward integrations return close to the initial state, which
  checks the time symmetry of the bridge composition.
- Earth-Moon convergence tests show the expected second-order trend when
  reducing `dt_outer`.
- Energy errors are observed as bounded oscillations rather than monotonic
  drift in the tested Earth-Moon setup.
- Long runs require backreaction from subsystem kicks to the host particle;
  one-way coupling produces unacceptable momentum drift.
- The Jupiter-Io-Europa-Ganymede example keeps the Laplace-angle diagnostic
  well behaved and serves as the reference multi-body moon-resonance example.
- The solar-system bridge example resolves both Earth-Moon and the inner
  Galilean system inside a Sun-plus-eight-planets main simulation, which is
  the current sandbox-scale validation configuration.

Long-term checks currently included in the repository focus on the solar-system
setup:

- `validation/solar_system_ias15/out_2000yr/` contains the standard 2000-year
  comparison figures.
- `validation/solar_system_ias15/out_20000yr/` contains the longer 20000-year
  stability run used to verify that the bridge shows bounded oscillatory error
  rather than secular breakup in this configuration.
- The bridge is not expected to match IAS15 phase-for-phase over very long
  spans; instead the important checks are bounded energy and angular-momentum
  error, stable Earth-Moon separation, and preserved resonant behavior.

IAS15 is used here as a reference trajectory for finite-time residuals and
diagnostics, not as the design target of the bridge method.  The comparison
workflow measures energy, angular momentum, Earth-Moon distance, Laplace
angle, and heliocentric phase/radius residuals for the bridged subsystems.

`dt_outer` controls how often main and subsystem simulations exchange
perturbations.  `dt_inner` controls the internal REBOUND integration of each
subsystem.  Once the subsystem is internally resolved, reducing `dt_inner`
further will not remove bridge coupling error; reducing `dt_outer` is then the
relevant accuracy knob.

### Solar-System Parameter Choice

The solar-system bridge setup has two resolved subsystems: Earth-Moon and
Jupiter-Io-Europa-Ganymede.  Parameter sweeps compare runtime, conserved
quantities, and the Jovian Laplace-angle residual against an IAS15 reference.

The current recommended balanced setting is:

```text
dt_outer_days    = 0.75
earth_moon_ratio = 3
jovian_ratio     = 12
dt_earth_moon    = 0.25 day
dt_jovian        = 0.0625 day
```

This is the setting used by `examples/solar_system_bridge.c`.  In the 20000-year
bridge-only sweep, it is slower than the most aggressive setting but reduces
the long-term Laplace-angle residual substantially.

Measured tradeoffs from the local 20000-year tests:

```text
setting         bridge time   max |Laplace residual|   note
1d / 5 / 15     185.95 s      16.48 deg                fastest tested stable run
0.75d / 3 / 10  195.49 s      10.20 deg                small runtime cost, better phase
0.75d / 3 / 12  209.97 s       7.28 deg                recommended balanced choice
0.75d / 3 / 15  233.68 s       5.02 deg                more conservative phase choice
```

Energy and angular-momentum errors remain at the same order of magnitude across
these runs.  The main accuracy tradeoff is long-term orbital phase, especially
the Jovian Laplace angle.  Larger outer steps such as `1.5d` are fast but showed
large Laplace residuals in the 5000-year screen, and non-integer-looking decimal
steps can require tolerant time synchronization because accumulated floating
point differences eventually reach a few `1e-12 yr`.

## Minimal Use

```c
#include "rebound.h"
#include "rebound_bridge.h"

struct reb_bridge* bridge = reb_bridge_create(main_sim, dt_outer);
reb_bridge_add_subsystem(bridge, host_index, sub_sim, dt_inner, 0);
reb_bridge_advance(bridge, 1.0);
reb_bridge_free(bridge);
```

## Examples

`examples/earth_moon_bridge.c` demonstrates bridge convergence for a resolved
Earth-Moon subsystem.

`examples/jovian_laplace_bridge.c` builds a resolved Jupiter-Io-Europa-Ganymede
subsystem from fixed JPL Horizons vectors for a stable Laplace-resonance
demonstration and prints the period ratios plus the Laplace angle diagnostic.

`examples/solar_system_bridge.c` builds a Sun plus eight-planet main simulation
with two resolved bridge subsystems: Earth-Moon and
Jupiter-Io-Europa-Ganymede. It is intended to validate multi-subsystem bridge
plumbing for solar-system-scale sandbox simulations. The example uses a larger
outer timestep for planetary motion and smaller inner timesteps for the resolved
moon subsystems.  The example uses the balanced sweep setting documented above:
`dt_outer = 0.75 day`, Earth-Moon ratio `3`, and Jovian ratio `12`.

## Validation Suite

`validation/solar_system_ias15/` contains a standalone bridge-vs-IAS15
comparison workflow for the solar-system setup.

- `compare_solar_system_ias15.c` runs both the bridge and a direct IAS15
  reference from the same initial conditions and writes a CSV of diagnostics.
- `plot_solar_system_ias15.py` renders energy, angular-momentum, distance, and
  phase/residual plots plus long-run stability and envelope plots from that CSV.
- `run.ps1` builds the validation target, runs the default 2000-year
  comparison, and generates the figures in
  `validation/solar_system_ias15/out_2000yr/`.
- `out_20000yr/` stores the longer 20000-year stability figures for the same
  solar-system configuration.
- `bench_whfast_bridge.c` and `plot_efficiency_sweep.py` provide the current
  runtime comparison workflow.
- `benchmark_figures/` stores the tracked bridge-vs-IAS15 and
  bridge-vs-WHFast speed comparison figures and CSV summaries.

## CMake Build

Configure with the default REBOUND 5 download:

```powershell
cmake -S . -B build
cmake --build build --config Release
```

Use a local REBOUND 5 source tree instead:

```powershell
cmake -S . -B build -DREBOUND_SRC_DIR="C:\path\to\rebound\src"
cmake --build build --config Release
```

Disable examples when building only the library:

```powershell
cmake -S . -B build -DREBOUND_BRIDGE_BUILD_EXAMPLES=OFF
cmake --build build --config Release
```

The CMake build creates static `rebound` and `rebound_bridge` targets. Link
downstream CMake code against `rebound_bridge::rebound_bridge` when consuming an
installed package export.

Example MSVC build from this directory:

```powershell
mkdir build
$env:REBOUND_SRC = "C:\path\to\rebound\src"
cl /nologo /Iinclude /I$env:REBOUND_SRC /Fo:build\ examples\jovian_laplace_bridge.c src\rebound_bridge.c $env:REBOUND_SRC\*.c /Fe:build\jovian_laplace_bridge.exe
build\jovian_laplace_bridge.exe
```

## License

`rebound_bridge` is licensed under GPL-3.0-or-later.  It builds against
REBOUND, which is also GPL-3.0-or-later.  The default CMake configuration
downloads REBOUND 5.0.0 from the upstream REBOUND repository and compiles it
into the local build.
