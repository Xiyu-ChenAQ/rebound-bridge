# Python API and compatibility contract

This reference describes the Python interface for **1.0.0b1**.
The API may change during the beta period. The 1.x compatibility policy below
starts with the final 1.0.0 release.
See the [quick start](../README.md#python-quick-start) for examples and installation.

## Public interface

Import public names from `rebound_bridge`: `Simulation`, `Body`, `Bridge`,
`Subsystem`, `BridgeError`, `CheckpointError` and `__version__`. Internal modules,
underscore-prefixed attributes, native pointers and serialization helpers are
implementation details. Direct construction of `Body` or `Subsystem`, and
subclassing these classes to override implementation methods, are outside this
contract; obtain handles through the factories below.

### Physical hierarchy: `Simulation` and `Body`

```python
Simulation(*, units=None, dt=None, dt_sub=None, integrator="whfast",
           method="lf4_2", coupling="monopole")
sim.add(name=None, **particle)                         # -> Body
body.add(name=None, **particle)                        # -> Body
sim.add_horizons(target, *, date=None, name=None, plane=None, m=None)  # -> Body
body.add_horizons(target, *, date=None, name=None, plane=None, m=None) # -> Body
sim.build()                                          # -> sim
sim.integrate(t, *, checkpoint=None, checkpoint_interval=None)
sim.to_rebound()                                      # -> rebound.Simulation
sim.save_to_file(filename)
Simulation.from_file(filename)                        # -> built Simulation
sim.hierarchy()                                      # -> str
sim.status()                                         # prints a summary
sim.close()
```

`add` accepts physical REBOUND particle keywords, including orbital elements.
A root uses global Cartesian coordinates; a child uses orbital elements or
Cartesian offsets relative to its parent. Choose the parent through `body.add`,
not a `primary` keyword. `name` sets particle identity and its REBOUND hash;
an explicitly supplied `hash` must agree. Names must be unique nonempty strings.
Root bodies cannot be added with orbital elements. `None` particle arguments
retain REBOUND's defaults: empty orbital values do not select orbital input,
and empty Cartesian values are zero parent-relative offsets. A positional
string in `add("name", ...)` is a local name, not a Horizons query.

`add_horizons` delegates the query to REBOUND's `Simulation.add`. Set units
before the first query and supply its explicit `date`. The initial reference
plane defaults to `"ecliptic"`; `"frame"` is also supported. Later queries
inherit the same date and plane, and conflicting explicit values are rejected.
`target` is a nonempty string. Dates accept `YYYY-MM-DD` with optional
`HH:MM[:SS]`, `JD` followed by a decimal Julian date, or a naive `datetime`.
`name` defaults to the target string; `m`, when provided, is in the simulation's mass
unit. Each imported state is already global, including when added through a
parent handle. The parent determines hierarchy membership, not a coordinate
offset. Query settings do not replace the simulation clock with a calendar or
Julian date. Use physical-body IDs when resolving satellites: e.g. Jupiter
`"599"`, not the Jupiter-system barycentre `"5"`. Existing global REBOUND
simulations can also be imported through `Bridge.from_simulation`.

| Member | Meaning and mutability |
| --- | --- |
| `sim.units` | REBOUND unit specification; set before adding bodies. Treat the returned mapping as a value, not an editing interface. |
| `sim.G` | Finite positive gravitational constant; set before build, after units. |
| `sim.t` | Current time; initial epoch may be assigned before build. |
| `sim.N` | Read-only number of physical bodies, excluding COM proxies. |
| `sim.horizons` | `None` before the first query; afterwards a read-only mapping with `date` and `plane`, preserved in tree checkpoints. |
| `sim.dt` | Outer coupling step; positive finite value required before build. |
| `sim.dt_sub` | Default internal step, or `None` if every resolved host specifies `dt`. |
| `sim.integrator` | Main integrator and inherited internal integrator; default `"whfast"`. |
| `sim.method`, `sim.coupling` | Coupling composition and interaction model, as defined below. |
| `sim.bodies`, `sim[name]` | Read-only name-to-`Body` mapping and name lookup. |
| `sim.children` | Tuple of root body handles. |
| `body.name`, `body.parent`, `body.children` | Read-only identity, parent handle (or `None`) and tuple of children. |
| `body.m` | Physical mass; editable before build. Finite and nonnegative; parents must have positive mass. |
| `body.dt` | Internal step override for a resolved host; `None` inherits `sim.dt_sub`. |
| `body.integrator` | Internal integrator override for a resolved host; otherwise inherits `sim.integrator`. |
| `body.effective_integrator` | Read-only integrator for the body's component; descendants share the resolved host's internal integrator, and bodies outside internal systems use the main integrator. A resolved host reports its internal integrator. |

All configuration setters freeze after a successful build. Changing `G` or mass
before build preserves existing Cartesian states; it does not recalculate
previously supplied orbital elements. Later additions use the updated values.
Only a root's direct child with descendants owns an internal simulation. All of
that child's descendants share it; the physical tree does not create recursive
Bridge time levels. Setting an internal step or integrator on another body is
rejected when building.

### Component interface: `Bridge` and `Subsystem`

```python
Bridge(simulation=None, *, dt, method="lf4_2", coupling="monopole")
bridge.add(particle=None, **kwargs)                    # -> None
bridge.add_subsystem(primary=None, *, dt, name=None,
                     integrator="whfast", **host)     # -> Subsystem
Bridge.from_simulation(simulation, *, subsystems, dt, dt_sub,
                       method="lf4_2", coupling="monopole",
                       integrator="whfast")           # -> unbuilt Bridge
subsystem.add(particle=None, **kwargs)                 # -> None
bridge.build()                                       # -> bridge
bridge.integrate(t, *, checkpoint=None, checkpoint_interval=None)
bridge.to_rebound()                                   # -> rebound.Simulation
bridge.save_to_file(filename)
Bridge.from_file(filename)                            # -> built Bridge
bridge.close()
```

The constructor copies an input `rebound.Simulation`; without one it creates
an empty WHFast main simulation with `main.dt = dt`. `from_simulation` takes a
mapping of subsystem names to disjoint nonempty lists of particle indices or
REBOUND hashes. The first member is the local host. `dt_sub` is a scalar or a
mapping from each subsystem name to its step; `integrator` selects the internal
integrator, while the main simulation retains the input integrator.

In `add_subsystem`, `primary` selects an existing main particle by index, hash
or particle handle. Alternatively, `**host` creates a new host using REBOUND
particle arguments; these two forms are mutually exclusive. A subsystem's first
particle is its host. Add members in the host's local Cartesian frame; the
orbital arguments default to that host, including when `primary=None` is
explicit. Orbital fields set to `None` do not activate this default for
Cartesian input. Subsystem Horizons name queries are rejected because they
return global states.

Assign each host a unique nonzero REBOUND hash before registering it if the
main particle array will be edited before build. The wrapper resolves that
identity again when building, so removing preceding particles or reordering
the main array does not redirect a subsystem to another host. Removing the
registered host or making its hash ambiguous is rejected. For a host registered
with a zero or duplicate hash, only an unchanged original slot/state with no
particle-count decrease can be verified; ambiguous edits are rejected. A
successful build updates `subsystem.host_index` to the resolved host index.

| Member | Meaning and mutability |
| --- | --- |
| `bridge.main` | Live, owned REBOUND main simulation; configurable before build within the supported model below. |
| `bridge.subsystems` | Read-only name-to-`Subsystem` mapping. |
| `bridge.t`, `bridge.dt` | Read-only current time and outer coupling step. |
| `bridge.method`, `bridge.coupling` | Selectors; the advanced interface permits changing them at completed outer boundaries. |
| `subsystem.name`, `subsystem.host_index` | Read-only identity and index of its main-system host/proxy. |
| `subsystem.simulation`, `subsystem.particles` | Live, owned internal REBOUND simulation and particle view. |
| `subsystem.dt`, `subsystem.integrator` | Internal settings; assign before build. |

Set the main component step through `bridge.main.dt` before build; it is
independent of `bridge.dt`. Set internal steps through `subsystem.dt`. A direct
assignment to `subsystem.simulation.dt` that conflicts with this setting is
rejected. Adaptive IAS15/BS timestep proposals may evolve during integration;
they are not changes to the configured nominal step.

## Lifecycle and ownership

`build()` validates and compiles the fixed hierarchy. Calling it again validates
the existing hierarchy and returns the same object. `integrate()` and
`save_to_file()` build automatically. `Simulation.to_rebound()` previews an
unbuilt tree without freezing it; **`Bridge.to_rebound()` builds and freezes the
advanced hierarchy**. Invalid target times and checkpoint cadences are checked
before an integration request builds the object.

Input simulations are copied and synchronized privately, leaving the caller's
states and caches unchanged. Building the advanced interface replaces component
simulations with copies in the appropriate frames: the main system contains
COM proxies and internal systems use local COM coordinates. Reacquire
`bridge.main` and `subsystem.simulation` after building; references obtained
earlier no longer address the live components.

After build, do not edit live component particles, masses, hashes, units, `G`,
integrators or fixed timesteps, and do not advance components independently.
Use `bridge.integrate()` to keep their times consistent. The wrapper catches
several invalid mutations, but arbitrary writes through raw REBOUND objects
are not a supported update mechanism. The advanced `method`/`coupling` setters
are explicit exceptions: changing them selects a different map/model for
subsequent steps and does not preserve a single unchanged integration scheme.

Particle counts, hashes and masses are checked exactly against the last accepted
state before public numerical operations and each native advance batch. A small
mass change may be accepted only when a successful native WHFast advance in
`barycentric` coordinates reconstructs its first particle's mass. This allowance
is bounded at both total-mass and first-particle-mass scales; it does not permit
user mass edits, even by one floating-point unit. Excessive reconstruction error
fails the live object and leaves the previous checkpoint available for restart.
For poorly conditioned mass ordering, use Jacobi coordinates or put the dominant
mass first. WHFast coordinates are fixed after build.

Every `to_rebound()` result is a **detached physical snapshot** in one global
frame, with each real body once and no COM proxies. Editing or integrating it
does not affect the hierarchy. It is not a substitute for a restart checkpoint:
it does not retain the Bridge hierarchy or component integration caches. Use
particle names/hashes for identity rather than relying on reconstructed order.

Both `Simulation` and `Bridge` support `with ... as ...` and `close()` for
releasing native resources. Do not reuse a closed object for numerical work.
Methods without a return annotation in the signatures above return `None`.
Restored objects are already built, with fixed membership and settings. Restore
tree checkpoints with `Simulation.from_file`; `Bridge.from_file` accepts both
tree and advanced checkpoints. An advanced checkpoint without tree metadata
cannot reconstruct a `Simulation` tree.

## Numerical and physical support

- `method="kdk"` selects second-order Bridge; `"lf4_2"` is the Python default.
- `coupling="monopole"` retains each subsystem's resolved interaction with the
  other subsystem's monopole. It omits interactions involving non-monopole
  moments of both systems. `"resolved"` includes all particle pairs between
  distinct subsystems. This model error differs from timestep error.
- The model is autonomous, unsoftened Newtonian point-mass gravity (`basic` or
  `compensated`) in fixed, weakly coupled hierarchies. There is no dynamic
  regrouping, collision handling, boundary handling or close-encounter adapter.
  Variational particles, custom `N_active` subsets, attached callbacks,
  REBOUNDx extras, particle `ap` data, user ODEs, MEGNO, live REBOUND displays
  and servers are unsupported. Native component Simulationarchive output is
  also unsupported; use Bridge checkpoints for restart. Massless
  members are allowed, but resolved hosts must have positive mass. Every
  nonempty WHFast or SABA component also requires its first particle to have
  positive mass; add the central massive body first or choose an appropriate
  alternative such as IAS15.
- Components must share units, `G` and epoch, and keep `exact_finish_time=1`.
  Steps must be positive and finite, with an integer outer/internal step ratio.
  Fixed-step component stages use equal substeps no larger than the configured
  component step. Choose steps by convergence checks for the intended observable.
- `integrate(t)` advances forward to an **absolute** target on the outer grid
  anchored at the initial epoch. It does not shorten a final step to hit an
  off-grid target. Adjacent grid boundaries must be distinguishable, and an
  advance must keep the outer step above the supported native time resolution
  at both the initial epoch and target. Fixed-step components must also retain
  resolvable half-step clock increments. Fixed stages are also rejected when
  accumulated clock roundoff cannot reliably preserve the planned number of
  substeps. A target
  may therefore be rejected even if adding a full outer step changes the time.
  Use a relative initial epoch and suitable units to keep the outer grid and
  fixed component steps above floating-point time resolution. Keep the
  original absolute epoch separately for output conversion.

Python IAS15 and BS components integrate each stage with a local time origin,
then restore the absolute component time. Their initial or adaptive timestep
may therefore be smaller than one floating-point unit of the absolute epoch.
This relies on the autonomous model above, whose forces do not depend explicitly
on absolute time; the outer-grid checks still apply.

Use these canonical integrator names in wrapper setters:

| Integrator | Supported configuration |
| --- | --- |
| `whfast` | Positive first-particle mass; `default` or `composition` kernel; correctors and synchronized/unsynchronized import are supported. `lazy` and `modifiedkick` are not adapted. |
| `ias15` | Ordinary adaptive N-body integration. |
| `bs` | Ordinary adaptive N-body integration, including its built-in N-body ODE; user ODEs or altered callbacks are excluded. |
| `leapfrog` | Ordinary fixed-step integration. |
| `saba` | Positive first-particle mass; `safe_mode=1`, `keep_unsynchronized=0`, uncorrected compositions only. |
| `eos` | `safe_mode=1`. |

TRACE, MERCURIUS, JANUS and other unadapted integrators are rejected. Aliases
that encode additional settings, such as `WHC` or `SABA(2)`, are rejected in
wrapper setters; configure supported advanced options before build. Choosing a
supported adaptive component does not make the overall map symplectic or
exactly time reversible. LF(4,2) retains its component stages even without
subsystems, so it need not reproduce a direct WHFast call's step sequence.

The current dependency is **REBOUND 4.6.0** and Python **3.10+**. Installation
from source requires a C compiler. CI has passed source builds, Python package
installation, paper analysis checks and short experiments on Windows, Linux
and macOS with Python 3.10 and 3.13. The full Python regression suite has also
been tested locally on Windows with Python 3.13.

## Checkpoint compatibility

`save_to_file(filename)` writes one restart snapshot at a completed outer
boundary. `integrate(..., checkpoint=filename, checkpoint_interval=interval)`
saves the initial state, each completed interval and the final state. The
interval is in simulation time and must be an integer multiple of the outer
step; its default is 1024 outer steps. Each successful write atomically replaces
the prior file. If native advancement or the following state update fails or
is interrupted, discard the live object and restore the last successfully
written checkpoint. After every batch, including the final one, the wrapper
validates the resulting state before saving or returning success. Non-finite
positions or velocities raise `ValueError`; the object is then marked failed
and cannot integrate, export or overwrite a checkpoint. Resolve the cause of
the failure before restarting the calculation.

Compatibility is determined by **format/schema and native state compatibility**,
not by the producer's `rebound-bridge` package version:

- The reader supports checkpoint `format_version=1`, hierarchy
  `bridge_schema=1` and, when a physical tree is present, tree `schema=1`.
  Missing or unsupported required schema versions are rejected before reading
  native simulation payloads. A checkpoint without a physical tree remains a
  valid advanced Bridge checkpoint.
- `bridge_version` records the producing package for provenance; it is not an
  exact-version gate. Existing 0.2.0 checkpoints with these schemas need no
  migration if their saved state passes validation. In particular, a checkpoint
  whose actual time differs from `epoch + completed_steps * dt` by more than
  `dt / 32` is rejected. Some older large-epoch calculations accumulated such
  clock drift. Read those files with their original producing version to
  inspect the saved state; the new reader does not modify the file or silently
  relabel its time. This is a maximum permitted drift, not a replacement for
  other state and epoch checks. Accepted smaller drift retains the saved
  physical state and actual component times when loading.
- The **exact REBOUND version**, recorded native ABI (including OS, architecture,
  byte order and relevant type sizes), payload checksums and hierarchy
  consistency must still match. Cross-platform checkpoint portability is not
  promised. Load only trusted files: REBOUND parses native binary payloads.
- Future readers must explicitly support a schema or provide a documented
  migration. An old reader is not required to read newer files; in particular,
  older 0.2.0 readers retain their original exact package-version check.

Format compatibility means the stored model and state can be restored; it does
not guarantee bitwise-identical trajectories across package releases, compilers
or changed call/output schedules. Native BS restart recreates internal work
arrays and likewise does not promise bitwise-identical continuation. Record
software versions and rerun accuracy checks when upgrading numerical software.

## Errors and 1.x version policy

`BridgeError` subclasses `RuntimeError` and reports native integration failures.
`CheckpointError` subclasses `ValueError` and reports invalid or incompatible
restart files. Invalid values and unsupported configurations normally raise
`ValueError`; wrong argument types may raise `TypeError`; lifecycle misuse
raises `RuntimeError`; allocation failure raises `MemoryError`. Unknown
attributes on tree `Simulation`/`Body` objects raise `AttributeError`.
Component integration errors preserve the upstream REBOUND failure diagnostic.
Catch exception classes rather than matching their messages.

Beginning with final 1.0.0, the documented imports, signatures, property
semantics, ownership and lifecycle above are the stable Python API for 1.x.
Patch releases may fix defects; minor releases may add backward-compatible
features. API removals, incompatible signature changes and changes to documented
defaults require a major release. Planned removals will be announced in release
notes and, where practical, emit `DeprecationWarning` for at least one minor
release before removal. Fixes that reject previously accepted invalid states
may ship without preserving that invalid behaviour and will be documented.

This promise does not freeze internal modules, binary layouts, exact floating-
point results, performance, exception text, `repr`/tree formatting, arbitrary
REBOUND attribute mutations or dependencies' own APIs. Changes to the supported
REBOUND pin or platform requirements must be documented; checkpoint migration
and continued support of earlier formats are stated separately. `__version__`
is the installed distribution's version string.
