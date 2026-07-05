# BCH Notes and Current Bridge Limitations

This note explains why the current `rebound_bridge` implementation should not
claim a fully BCH-derived high-order force-gradient method.

The short version is:

- BCH analysis is useful for understanding ideal Hamiltonian splittings.
- The current Bridge code uses host-relative subsystem variables.
- Those variables are not written as one explicit global canonical coordinate
  system.
- Therefore the standard Cartesian BCH force-gradient coefficients cannot be
  applied blindly to the current implementation.

This does not mean BCH theory is wrong.  It means the implemented variables do
not yet satisfy the assumptions needed for the standard textbook derivation.

## What BCH Assumes

The usual BCH derivation assumes that the split operators are Hamiltonian flows
on the same canonical phase space.  In practice that means:

- both parts of the split act on the same coordinates and momenta;
- the coordinates are canonical;
- the drift and kick maps are exact flows of the split Hamiltonians;
- the kinetic energy has the usual diagonal Cartesian mass metric.

Those assumptions are clean in a single global Cartesian system.  They are not
automatic in a Bridge implementation that stores each subsystem in a
host-relative frame.

## What the Current Bridge Does

The current implementation coordinates separate REBOUND simulations:

- a main simulation containing the Sun, planets, and host particles;
- subsystem simulations containing resolved members in host-relative frames;
- cross-kicks between the main simulation and the resolved subsystem members.

For a subsystem member, the global state is reconstructed from the host state
plus the member's relative state.

This is a practical and useful representation, but it is not the same thing as
declaring one global canonical Hamiltonian split.  In relative variables, the
kinetic energy generally has a non-diagonal mass metric.  The simple Cartesian
force-gradient expression assumes the diagonal case.

That is the main limitation: the force-gradient vector computed by the code can
be internally consistent with the programmed cross-kick, while still not being
the canonical BCH commutator required by the standard high-order theorem.

## What Was Tested

Several higher-order ideas were tested:

- the plain second-order Bridge KDK step;
- the EOS LF4_2 composition;
- an experimental forward force-gradient fourth-order method;
- a finite-difference check of the force-gradient acceleration.

The finite-difference check confirms that the analytic force-gradient vector
matches the directional derivative of the current C cross-acceleration map.
For one Pluto-Charon diagnostic run, the analytic and finite-difference values
agreed to within about `3.5e-12`, while a global sign error would have produced
a difference around `3.8e-10`.

![Finite-difference force-gradient diagnostic](bch_force_gradient_fd_check.png)

This rules out a simple sign mistake, but it does not prove that the current
map is the canonical BCH force-gradient operator.

## Numerical Evidence

In 100-year Pluto-Charon tests with `dt_outer = 2 days`, the approximate maximum
phase-arc errors were:

- order2: about `15.2 km`;
- EOS LF4_2: about `0.30 km`;
- forward FG4: about `1.23 km`.

With `dt_outer = 4 days`, EOS LF4_2 and the forward force-gradient method no
longer showed a clean fourth-order-style advantage.

![Pluto-Charon method comparison](bch_pluto_charon_phase_methods.png)

The force-gradient scale was also swept over a deliberately large range.  The
phase error barely changed, which suggests that this channel is not the main
control knob for the remaining fast-binary phase drift.

![Forward FG4 force-gradient scale sweep](bch_force_gradient_scale_sweep.png)

## Current Interpretation

The correct interpretation is conservative:

- The current Bridge implementation is useful and stable in the tested cases.
- EOS LF4_2 gives a real improvement for fast-binary phase behavior.
- The experimental force-gradient derivative is internally consistent with the
  current cross-kick map.
- The current code should not claim that forward FG4 solves the fast-binary
  phase problem.
- The current code should not claim that the standard BCH coefficients apply
  directly to the host-relative representation.

This is why the force-gradient code remains in `_experimental.c` files.

## What Would Be Needed

A rigorous BCH-based high-order Bridge would need one of these routes:

1. Rewrite the Bridge in explicit canonical variables, such as a Jacobi-like or
   democratic-heliocentric-like coordinate system, and derive the correct mass
   metric.
2. Build a global Cartesian prototype that expands all subsystem members into
   one phase space, applies the high-order kernel there, and then studies how
   to project back safely.
3. Use the current EOS LF4_2 method as the practical default and treat any
   remaining fast-binary phase calibration as empirical rather than BCH-proven.

For this repository, the recommended method is currently
`REB_BRIDGE_METHOD_EOS_LF4_2`.

The open research question is whether a resolved-subsystem Bridge can be
formulated in explicit canonical variables so that BCH, force-gradient, or
processor theory applies directly to fast hierarchical binaries.
