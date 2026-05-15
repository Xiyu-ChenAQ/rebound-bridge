# rebound_bridge

Standalone symplectic bridge scheduler for REBOUND.

It only depends on the REBOUND C API and can be moved into its own repository or built as a
separate static/shared library later.

## Layout

```text
rebound_bridge/
  include/rebound_bridge.h
  src/rebound_bridge.c
  examples/earth_moon_bridge.c
  examples/jovian_laplace_bridge.c
  examples/solar_system_bridge.c
```

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

## Validation Notes

`rebound_bridge` is a scheduler around REBOUND simulations, not a standalone
replacement integrator.  A typical configuration uses WHFast inside the main
simulation and inside each subsystem, while the bridge supplies the symmetric
KDK coupling between subsystem barycenters and resolved subsystem bodies.

The current validation focus is long-term structural behavior:

- Forward-then-backward integrations return close to the initial state, which
  checks the time symmetry of the bridge composition.
- Earth-Moon convergence tests show the expected second-order trend when
  reducing `dt_outer`.
- Energy errors are observed as bounded oscillations rather than monotonic
  drift in the tested Earth-Moon setup.
- Long runs require backreaction from subsystem kicks to the host particle;
  one-way coupling produces unacceptable momentum drift.
- IAS15 runs are used as reference trajectories for measuring residuals and
  phase error over finite spans, not as the design target of the bridge method.

`dt_outer` controls how often main and subsystem simulations exchange
perturbations.  `dt_inner` controls the internal REBOUND integration of each
subsystem.  Once the subsystem is internally resolved, reducing `dt_inner`
further will not remove bridge coupling error; reducing `dt_outer` is then the
relevant accuracy knob.

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
moon subsystems.

## Validation Suite

`validation/solar_system_ias15/` contains a standalone bridge-vs-IAS15
comparison workflow for the solar-system setup. The default workflow now runs
a 2000-year long-term stability check.

- `compare_solar_system_ias15.c` runs both the bridge and a direct IAS15
  reference from the same initial conditions and writes a CSV of diagnostics.
- `plot_solar_system_ias15.py` renders energy, angular-momentum, distance, and
  phase/residual plots plus long-run stability and envelope plots from that CSV.
- `run.ps1` builds the validation target, runs the 2000-year comparison, and
  generates the figures in `validation/solar_system_ias15/out/`.

## CMake Build

Configure with a REBOUND C source tree:

```powershell
cmake -S . -B build -DREBOUND_SRC_DIR="C:\path\to\rebound\src"
cmake --build build --config Release
```

If `_deps/rebound-upstream/src` exists, CMake uses it as the default
`REBOUND_SRC_DIR`.

Disable examples when building only the library:

```powershell
cmake -S . -B build -DREBOUND_SRC_DIR="C:\path\to\rebound\src" -DREBOUND_BRIDGE_BUILD_EXAMPLES=OFF
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
