# rebound_bridge

`rebound_bridge` couples REBOUND simulations that evolve on different
timescales. It is designed for fixed, weakly coupled hierarchies, such as
planets with compact satellite systems: each internal system takes short steps
while the main simulation follows its orbit with longer steps. The coupling
runs in C99 and is available through C and Python interfaces.

This is **version 1.0.0b1**, the first beta towards the 1.0 release.
Release notes are available on [GitHub Releases](https://github.com/Xiyu-ChenAQ/rebound-bridge/releases).

## Python quick start

From the repository root, install with:

```sh
python -m pip install .
```

This requires Python 3.10+ and a C compiler. The Python package pins REBOUND
4.6.0. Then create a star, a planet and two moons:

```python
import rebound_bridge

sim = rebound_bridge.Simulation()
sim.units = ("yr", "AU", "Msun")
sim.dt = 0.002       # Outer coupling step, in years.
sim.dt_sub = 0.0002  # Internal step for the planet and its moons.

star = sim.add(name="star", m=1.0)
planet = star.add(name="planet", m=1e-3, a=1.0)
planet.add(name="moon1", m=1e-7, a=0.01)
planet.add(name="moon2", m=2e-7, a=0.02, e=0.01, f=1.0)

print(sim.hierarchy())
sim.integrate(1.0)
physical = sim.to_rebound()
print(physical.t, physical.N)
```

The hierarchy follows the physical system: `sim.add(...)` adds a root body,
and `parent.add(...)` adds a body orbiting that parent. Bridge constructs the
centre-of-mass proxies and internal simulations automatically. The default
coupling method is LF(4,2), with WHFast for the component integrations. The
numbers above illustrate usage; choose steps through convergence checks for
your own system.

### Add bodies

Use familiar REBOUND particle arguments: `m`, `a`, `e`, `inc`, `Omega`, `omega`
and `f`. Orbital elements are relative to the parent receiving `.add(...)`,
so no `primary=` argument is needed. Angles are in radians, measured in the
common coordinate axes. Alternatively, provide `x`, `y`, `z`, `vx`, `vy`, `vz`;
for a child these are offsets from its parent, and for a root they are global
coordinates. Each `name` also becomes the particle's REBOUND hash.

To add another planet with its own satellite, insert this before the first
`sim.integrate(...)` call:

```python
other = star.add(name="other_planet", m=3e-4, a=2.0, e=0.02)
other.add(name="other_moon", m=1e-8, a=0.008)
```

Set units before adding bodies, and finish adding bodies and choosing settings
before calling `build()`, `integrate()` or `save_to_file()`, which freeze the
configuration. Each planet with moons shares one internal simulation with all
its descendants; planets without descendants remain in the main system.
Deeper nesting does not introduce additional Bridge timestep levels.

For ephemeris initial states, `sim.add_horizons(...)` and
`parent.add_horizons(...)` query JPL Horizons through REBOUND. Set units and
provide a date on the first query; later bodies reuse that date. See the
[API reference](python/API.md#physical-hierarchy-simulation-and-body) for query
arguments and coordinate conventions.

### Choose timesteps and integrators

These optional settings also go before the first integration:

```python
sim.method = "lf4_2"       # Or "kdk" for second-order Bridge.
sim.integrator = "whfast"  # Main integrator and default for internal systems.
sim.dt = 0.002            # Outer coupling step.
sim.dt_sub = 0.0002       # Default internal step.
planet.dt = 0.0001        # Override for the planet–moon internal system.
planet.integrator = "ias15"
```

The outer step must be an integer multiple of each configured internal step.
IAS15 adapts its steps within each component evolution. `planet.integrator`
selects the integrator for the planet and its moons; the planet system's orbit
around the star is advanced by the main integrator.
Use `body.effective_integrator` to inspect the integrator of its component,
including for a moon that shares its planet's internal simulation.

`sim.integrate(t)` advances to an **absolute target time**, in the chosen time
unit. The target must lie on the outer-step grid: after `integrate(1.0)`,
`integrate(2.0)` advances by one more year in this example.

### Read positions and orbits

`sim.to_rebound()` returns an ordinary `rebound.Simulation` containing each
physical body once, in a common inertial frame, without proxy particles. Use
REBOUND's particle and orbit syntax on this snapshot:

```python
physical = sim.to_rebound()
ps = physical.particles
print(ps["moon1"].x, ps["moon1"].y, ps["moon1"].z)
orbit = ps["moon1"].orbit(primary=ps["planet"])
print(orbit.a, orbit.e, orbit.inc)
```

The snapshot is independent: editing it does not change the Bridge simulation.
`sim.N` gives the physical body count, and `sim.status()` prints a brief summary.

### Save and resume

Save the complete hierarchy, settings and integration state:

```python
sim.save_to_file("run.rbb")
sim.close()
```

Resume in a new Python session without rebuilding the hierarchy:

```python
import rebound_bridge

with rebound_bridge.Simulation.from_file("run.rbb") as sim:
    sim.integrate(2.0, checkpoint="run.rbb", checkpoint_interval=0.1)
```

Here the checkpoint is updated every 0.1 years of simulation time and at the
final target. Its interval must be an integer multiple of `sim.dt`. After an
interruption, load the last completed checkpoint with `from_file()` again.

## Methods and limitations

The default Python method, `lf4_2`, uses the weak residual interaction to
reduce outer splitting error. It can show near-fourth-order convergence in
weakly coupled systems, but is not unconditionally fourth order. The `kdk`
method provides the standard second-order Bridge composition. Both use
REBOUND integrators to advance the components; check convergence in both
outer and internal steps. Once the internal motion is resolved, decreasing
`dt_sub` alone does not remove the outer splitting error.

With multiple subsystems, the default coupling resolves each subsystem
against the other's centre-of-mass proxy. It omits interactions involving
non-monopole structure in both subsystems; the leading omitted term is generally
quadrupole–quadrupole. This model difference remains as timesteps decrease. Set
`sim.coupling = "resolved"` before integration to include all member pairs
between distinct subsystems. A single resolved subsystem has no such
inter-subsystem approximation.

WHFast and IAS15 are the usual component choices; other supported integrators
and their restrictions are listed in the [API reference](python/API.md#numerical-and-physical-support).
TRACE, MERCURIUS and JANUS are not yet adapted by the Python layer. Subsystem
membership stays fixed: encounters that disrupt the hierarchy are outside the
intended regime. For ordinary planetary systems without fast internal
subsystems, direct REBOUND integration is usually the simpler choice.

## C interface

From a repository checkout or extracted source archive, build the C library
with CMake 3.20+ and a C99 compiler. Both include the REBOUND C sources.

```sh
cmake -S . -B build
cmake --build build --config Release
```

The public C interface is declared in
[rebound_bridge.h](src/include/rebound_bridge.h). Link CMake targets to
`rebound_bridge::rebound_bridge`.

## Testing

The public test suite contains the paper experiments and their supporting
checks. From the source root:

```sh
python -m pip install -r tests/paper/requirements.txt
python tests/paper/run_all.py --unit-tests-only
python tests/paper/run_all.py --mode smoke --output-dir build/paper-smoke
```

The first command after installation checks the experiment and plotting scripts.
The smoke run builds the C experiment programs, runs their registered checks,
and exercises each experiment with shortened settings. Use a fresh output
directory for each run. These checks verify execution and output; reproducing
the seven manuscript figures requires the full settings in the
[paper guide](tests/paper/README.md).

## Platforms

Version 1.0.0b1 has been tested locally on Windows with Python 3.13
and REBOUND 4.6.0. CI has also passed on Windows, Linux and macOS with Python
3.10 and 3.13, covering C builds, Python package installation, paper analysis
checks and short experiments. These checks do not replace convergence tests
for your own system.

## Documentation

- [Python quick start](#python-quick-start): a first simulation, output and restart.
- [API reference](python/API.md): supported configurations, ownership and checkpoint compatibility.
- [C API header](src/include/rebound_bridge.h): declarations for the low-level interface.
- [Testing](#testing) and [paper experiments](tests/paper/README.md): validation and numerical reproduction.

The source distribution includes the paper experiments, their initial
conditions and plotting scripts. Use the source archive for reproduction;
the wheel installed by `pip` contains the library. The paper guide maps all
seven figures to individual experiments and provides a short smoke run.

The lower-level Python `Bridge` interface accepts existing REBOUND simulations;
see the [component interface](python/API.md#component-interface-bridge-and-subsystem).

## Feedback

This is a beta release of `rebound_bridge`. If you run into a
bug, find the documentation confusing, or wish it could do something it doesn't
yet, please [open an issue](https://github.com/Xiyu-ChenAQ/rebound-bridge/issues).
Questions and suggestions are welcome too. I'd be glad to hear what you're
trying to simulate.

## Citation

Software citation metadata is provided in [CITATION.cff](CITATION.cff).
Please also cite the associated methods paper once it is available.

## License

`rebound_bridge` is licensed under [GPL-3.0-or-later](LICENSE). The vendored
REBOUND code uses the same license.
