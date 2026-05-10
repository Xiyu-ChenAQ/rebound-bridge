# rebound_bridge

Standalone symplectic bridge scheduler for REBOUND.

This directory is intentionally independent from `cs/`.  It only depends on
the REBOUND C API and can be moved into its own repository or built as a
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

Example MSVC build from this directory:

```powershell
mkdir build
$env:REBOUND_SRC = "C:\path\to\rebound\src"
cl /nologo /Iinclude /I$env:REBOUND_SRC /Fo:build\ examples\jovian_laplace_bridge.c src\rebound_bridge.c $env:REBOUND_SRC\*.c /Fe:build\jovian_laplace_bridge.exe
build\jovian_laplace_bridge.exe
```
