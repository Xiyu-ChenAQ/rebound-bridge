# Pluto-Charon Binary Validation Design

## Purpose

Validate `rebound_bridge` on a true hierarchical binary system, not only a planet-satellite system. The Earth-Moon barycenter remains inside Earth, while the Pluto-Charon barycenter is outside Pluto. A stable Pluto-Charon result would show that EOS Bridge can handle resolved binary subsystems whose main-simulation host particle represents the subsystem center of mass, not a physical primary body.

## Scope

Implement a new validation area:

```text
validation/pluto_charon_binary/
```

The validation will have two stages:

1. Minimal case: `Sun + Pluto-Charon`.
2. Extended case: `Sun + Jupiter + Saturn + Uranus + Neptune + Pluto-Charon`.

The minimal case isolates Bridge splitting and coordinate handling. The extended case adds realistic outer-solar-system perturbations after the minimal case is working.

## Coordinate Model

The main simulation must contain a single host particle for the Pluto-Charon binary. That host is the Pluto-Charon center of mass and has mass `m_pluto + m_charon`.

The subsystem contains two physical particles, Pluto and Charon, expressed in the Pluto-Charon COM frame:

```text
sum(m_i * r_i) = 0
sum(m_i * v_i) = 0
```

Bridge snapshots reconstruct global inertial Pluto and Charon states as:

```text
r_pluto_global  = r_pc_com_host + r_pluto_sub
r_charon_global = r_pc_com_host + r_charon_sub
v_pluto_global  = v_pc_com_host + v_pluto_sub
v_charon_global = v_pc_com_host + v_charon_sub
```

The IAS15 reference simulation contains physical Sun, Pluto, and Charon particles in the same global inertial frame. For the extended case, the IAS15 reference also includes the outer planets as physical particles.

Pluto must never be used as the Bridge host for this validation. Doing that would test a planet-satellite approximation, not a true binary COM host.

## Initial State Invariants

The validation executable will support a zero-time initial-state check before any integration. It will compare reconstructed Bridge global states against the IAS15 physical-particle states.

Required initial-state diagnostics:

- `max_abs_mass_diff`
- `max_dr_au`
- `max_dv_au_per_yr`
- Pluto row: mass, position, velocity, `dr`, `dv`
- Charon row: mass, position, velocity, `dr`, `dv`
- Pluto-Charon COM row: reconstructed Bridge COM state, IAS15 mass-weighted COM state, `dr`, `dv`

The expected initial-state error is machine-level, approximately `1e-14` to `1e-13` in AU and AU/year depending on coordinate magnitude.

## Runtime Diagnostics

Each output sample will record both external COM behavior and internal binary behavior.

External binary orbit diagnostics:

- Pluto-Charon COM absolute position error versus IAS15.
- Pluto-Charon COM velocity error versus IAS15.
- Pluto-Charon COM orbital elements around the Sun when the Sun is present.

Internal binary diagnostics:

- Pluto absolute position error versus IAS15.
- Charon absolute position error versus IAS15.
- Pluto-Charon relative position error.
- Pluto-Charon distance difference.
- Pluto-Charon relative semimajor axis difference.
- Pluto-Charon relative eccentricity difference.

Conservation diagnostics:

- Bridge relative energy error.
- IAS15 relative energy error.
- Bridge relative angular momentum error.
- IAS15 relative angular momentum error.
- Bridge relative total momentum error.
- IAS15 relative total momentum error.

Timing diagnostics:

- Bridge step time and cumulative time.
- IAS15 step time and cumulative time.
- Overall speedup versus IAS15.
- Optional fair-speed comparison against pure WHFast, using the smallest timestep used inside Bridge.

## Error Attribution

The COM error is a first-class diagnostic, not a secondary plot.

If COM error remains small while Pluto-Charon relative orbit error grows, the dominant issue is likely subsystem-internal evolution or operator splitting inside the resolved binary.

If COM error grows significantly, the dominant issue is more likely the external coupling: cross-kick, host-proxy handling, or the main-system interaction model.

This separation is required because absolute Pluto and Charon errors combine two effects:

```text
absolute body error = binary COM error + internal relative-orbit error
```

Without COM diagnostics, the validation cannot cleanly distinguish external orbital drift from internal binary phase drift.

## Method Sweep

The first sweep will compare:

- `order2`
- `eos-lf4-2`

The primary timestep sweep will vary the outer Bridge timestep. The inner binary timestep ratio will be fixed for the first pass, then varied only if the COM error is small but internal orbit error is large.

The expected scientific question is:

```text
Can EOS LF4_2 preserve a strong external timestep advantage while keeping both
Pluto-Charon COM motion and Pluto-Charon relative orbit accurate?
```

## Outputs

The validation should produce:

```text
validation/pluto_charon_binary/out/
validation/pluto_charon_binary/benchmark_figures/
```

Core CSV files:

- Per-sample comparison CSV for the minimal case.
- Per-sample comparison CSV for the extended outer-planets case.
- Method and timestep sweep summary CSV.
- Initial-state CSV.

Core plots:

- COM error versus time.
- Pluto and Charon absolute error versus time.
- Pluto-Charon relative error versus time.
- Pluto-Charon orbital element residuals.
- Conservation diagnostics.
- Speedup versus error tradeoff.
- Timestep convergence plot.

## Acceptance Criteria

The validation is considered implemented when:

- The zero-time Bridge reconstructed states match IAS15 physical states at machine-level tolerance.
- Minimal and extended cases both run with `order2` and `eos-lf4-2`.
- The output clearly separates COM error from Pluto-Charon relative orbit error.
- The sweep summary can identify whether failures come from COM drift or internal binary drift.
- A fair pure-WHFast speed baseline can be produced using Bridge's smallest internal timestep.

## Implementation Notes

The new validation should be independent of `validation/solar_system_ias15` at the executable level. Shared helper code may be extracted later only if duplication becomes a maintenance problem. The first implementation should favor explicit coordinate handling over premature abstraction.

All initial conditions should be frozen in source or checked-in data files. The validation must not depend on runtime network access.
