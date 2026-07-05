# rebound_bridge

`rebound_bridge` is a small C library for coupling multiple REBOUND
simulations with an operator-splitting Bridge scheme.

The library vendors the REBOUND C source tree under `rebound/src`; a normal
CMake build does not download or locate an external REBOUND installation.

## Model

The Bridge owns no particles by itself.  It coordinates:

- one main REBOUND simulation containing slow bodies and one host particle per
  resolved subsystem;
- one or more subsystem REBOUND simulations containing resolved bodies in each
  subsystem frame;
- symmetric cross-kicks between the main system and resolved subsystems.

For a fixed-timestep REBOUND integrator such as WHFast:

- `dt_outer` is the main simulation timestep and Bridge coupling interval;
- `dt_inner` is the subsystem timestep;
- `dt_outer / dt_inner` must currently be an integer ratio.

## Public Headers

The only installed public header is `include/rebound_bridge.h`.  It exposes
the stable Bridge scheduler API.

Example-only helper declarations live in `examples/rebound_bridge_examples.h`
and are not installed as public API.

## Source Layout

The implementation follows the short naming style used by REBOUND's `src/`
tree.  Files ending in `_experimental` contain higher-order or force-gradient
compositions that are still under investigation.

- `src/bridge.c` — lifecycle, configuration, and subsystem registration.
- `src/bridge_kick.c` — plain Bridge cross-kick kernel.
- `src/bridge_advance.c` — order-2 KDK advance and stable dispatcher.
- `src/bridge_eos.c` — EOS LF4_2 advance method.
- `src/bridge_fg4_experimental.c` — forward force-gradient fourth-order advance.
- `src/bridge_force_gradient_experimental.c` — force-gradient cross-kick kernel
  and finite-difference diagnostics.
- `src/bridge_processor.c` — commutator processor.
- `src/bridge_internal.h` — private to the library implementation.
- `examples/bridge_examples.c` — convenience constructors used only by the
  bundled examples; not compiled into the installed library.

## Build

Default CMake build, which builds the vendored REBOUND static library and the
Bridge static library:

```sh
cmake -S . -B build
cmake --build build --config Release
```

The build creates two static targets:

- `rebound_bridge::rebound`, built from vendored `rebound/src/*.c`;
- `rebound_bridge::rebound_bridge`, built from `src/bridge*.c` and linked
  against `rebound_bridge::rebound`.

Optional examples, validation tools, and tests:

```sh
cmake -S . -B build -DREBOUND_BRIDGE_BUILD_EXAMPLES=ON -DREBOUND_BRIDGE_BUILD_VALIDATION=ON -DREBOUND_BRIDGE_BUILD_TESTS=ON
```

Downstream CMake projects should link against:

```cmake
target_link_libraries(my_target PRIVATE rebound_bridge::rebound_bridge)
```

Manual MSVC build from this directory on Windows:

```powershell
mkdir build
cl /nologo /Iinclude /Iexamples /Irebound\src ^
  examples\solar_system_bridge.c examples\bridge_examples.c src\*.c rebound\src\*.c ^
  /Fo:build\ /Fe:build\solar_system_bridge.exe
```

On Unix-like systems, prefer CMake.  If building manually, use the same include
and source sets with your C compiler and link the math library when needed.

## Platform Support

The core library is intended to be platform-neutral C99.  The library sources
under `src/` do not include `windows.h` and do not require Windows APIs.

Platform-specific code is limited to guarded support paths:

- `include/rebound_bridge.h` uses `__declspec(...)` only when building a shared
  library on Windows.
- `validation/validation_timer.h` uses `QueryPerformanceCounter` on Windows and
  `clock_gettime` or `timespec_get` elsewhere.

## Minimal Use

```c
#include "rebound.h"
#include "rebound_bridge.h"

struct reb_bridge* bridge = reb_bridge_create(main_sim, dt_outer);
reb_bridge_set_method(bridge, REB_BRIDGE_METHOD_EOS_LF4_2);
reb_bridge_add_subsystem(bridge, host_index, sub_sim, dt_inner, 0);
reb_bridge_advance(bridge, duration);
reb_bridge_free(bridge);
```

The helper functions `reb_bridge_set_integrator_whfast()` and
`reb_bridge_set_integrator_ias15()` configure the corresponding REBOUND
integrator on an existing `struct reb_simulation*`.

## Bridge Methods

- `REB_BRIDGE_METHOD_ORDER2` — symmetric KDK Bridge (2nd order).
- `REB_BRIDGE_METHOD_EOS_LF4_2` — embedded-operator-splitting LF4_2 (4th order).
- `REB_BRIDGE_METHOD_EOS_LF8_6_4` — *experimental*; higher-order EOS composition
  using LF8_6_4 coefficients.  Includes negative substeps; treat as a validation
  target rather than a default recommendation.
- `REB_BRIDGE_METHOD_FORWARD_FG4` — *experimental*; forward-only 4th-order
  force-gradient Bridge kernel.

Advanced cross-kick and corrector controls (`REB_BRIDGE_CROSS_KICK_FORCE_GRADIENT`,
`REB_BRIDGE_CORRECTOR_COMMUTATOR_AB`, `reb_bridge_force_gradient_finite_difference_check`)
are also declared in `rebound_bridge.h`.  Their coefficients and generality are
still under investigation; treat them as experimental.

## Examples

Bundled examples are enabled by `REBOUND_BRIDGE_BUILD_EXAMPLES=ON`:

- `examples/earth_moon_bridge.c` uses the example convenience constructors for
  a simple Earth-Moon Bridge demonstration.
- `examples/jovian_laplace_bridge.c` builds a resolved
  Jupiter-Io-Europa-Ganymede subsystem and prints Laplace-resonance diagnostics.
- `examples/solar_system_bridge.c` builds a Sun plus eight-planet main
  simulation with Earth-Moon and inner Galilean resolved subsystems.

Example-only declarations are in `examples/rebound_bridge_examples.h` so the
public `include/` directory stays focused on library API.

## Validation

Validation tools are enabled by `REBOUND_BRIDGE_BUILD_VALIDATION=ON`.

Solar-system validation:

- `validation/solar_system_ias15/compare_solar_system_ias15.c` compares Bridge
  against a direct IAS15 reference and writes diagnostics to CSV.
- `validation/solar_system_ias15/bench_whfast_bridge.c` benchmarks Bridge and
  WHFast configurations.

Earth-Moon / Galilean validation:

- `validation/earth_moon_whfast_bridge/compare_earth_moon_whfast_bridge.c`
  compares Bridge, WHFast, and IAS15 on a solar system with two resolved
  subsystems (Earth-Moon and inner Galilean moons).

Pluto-Charon validation:

- `validation/pluto_charon_binary/compare_pluto_charon.c` validates Bridge on a
  true binary subsystem.
- `validation/pluto_charon_binary/bench_pluto_charon.c` benchmarks the same
  setup.

Run the test suite:

```sh
cmake -S . -B build-pluto-charon -DREBOUND_BRIDGE_BUILD_TESTS=ON
cmake --build build-pluto-charon --config Release
ctest --test-dir build-pluto-charon -C Release --output-on-failure
```

The CTest suite includes a forward-backward time-reversibility check for the
public `reb_bridge_step()` KDK step.

## Notes

`rebound_bridge` is a scheduler around REBOUND simulations, not a standalone
replacement integrator.  The usual workflow is to use WHFast or IAS15 inside
the main and subsystem simulations while the Bridge supplies the cross-system
operator splitting.

The most important accuracy knob is generally `dt_outer`, because it controls
how often the main system and subsystems exchange perturbations.  Once a
subsystem is internally resolved, reducing `dt_inner` further will not remove
Bridge splitting error.

## License

`rebound_bridge` is licensed under GPL-3.0-or-later.  It vendors and builds
against REBOUND, which is also GPL-3.0-or-later.
