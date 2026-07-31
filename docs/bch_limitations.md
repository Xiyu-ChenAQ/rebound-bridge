# BCH Notes and Current Bridge Limitations

This note explains why the current `rebound_bridge` implementation should not
claim a fully BCH-derived high-order force-gradient method.

The short version is:

- BCH analysis is useful for understanding ideal Hamiltonian splittings.
- The current Bridge code uses host-relative subsystem variables.
- Those variables define a constrained symplectic representation for the
  ordinary Bridge split, but they are not independent particle-wise canonical
  pairs.
- Therefore the standard Cartesian BCH force-gradient coefficients cannot be
  applied blindly to the current implementation.

This limitation concerns the experimental force-gradient operator.  It does
not invalidate the position-only residual kick or the EOS LF4_2 composition.

## What BCH Assumes

The usual BCH derivation assumes that the split operators are Hamiltonian flows
on the same canonical phase space.  In practice that means:

- both parts of the split act on the same coordinates and momenta;
- the coordinates are canonical;
- the drift and kick maps are exact flows of the split Hamiltonians.

The commonly used Cartesian force-gradient formula makes an additional
assumption: the kinetic energy has the usual diagonal mass metric.  That
simplification is automatic in inertial Cartesian coordinates but not in a
redundant, constrained set of internal coordinates.

## What the Current Bridge Does

The current implementation coordinates separate REBOUND simulations:

- a main simulation containing the Sun, planets, and host particles;
- subsystem simulations containing resolved members in host-relative frames;
- cross-kicks between the main simulation and the resolved subsystem members.

For a subsystem member, the global state is reconstructed from the host state
plus the member's relative state.

Let the host carry the subsystem centre-of-mass variables, and let the internal
positions and momenta satisfy their mass-weighted constraints.  The canonical
one-form then separates into centre-of-mass and internal terms.  Thus the
constraint surface supplies the symplectic phase space needed for the ordinary
Hamiltonian split.  The residual cross Hamiltonian depends only on positions,
so its flow is the velocity kick implemented by the code.  Registration checks
the constraints; the integrator does not repeatedly project the subsystem back
onto them.

This result is sufficient for KDK and LF4_2.  It does not by itself identify the
correct force-gradient commutator in the redundant internal variables.  After
eliminating the constraints, the independent internal coordinates need not
retain the diagonal Cartesian mass metric assumed by the simple formula.

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
- The position-only residual is an exact kick on the constrained phase space.
- EOS LF4_2 is a valid outer composition; it is symplectic when its component
  maps are symplectic, and its generalized-order scaling requires component
  integration errors to remain below the outer error.
- The experimental force-gradient derivative is internally consistent with the
  current cross-kick map.
- The current code should not claim that forward FG4 solves the fast-binary
  phase problem.
- The current code should not claim that the standard BCH coefficients apply
  directly to the host-relative representation.

This is why the force-gradient code remains in `_experimental.c` files.

## What Would Be Needed

A rigorous BCH-based high-order Bridge would need one of these routes:

1. Derive the force-gradient commutator directly on the constrained phase
   space, including the projected inverse mass metric.
2. Introduce independent hierarchical Jacobi-like coordinates and derive the
   same operator there.
3. Build a global Cartesian prototype, where the standard mass metric applies,
   and compare its force-gradient map with the constrained implementation.

For this repository, the recommended method is currently
`REB_BRIDGE_METHOD_EOS_LF4_2`.

The open research question is whether the experimental force-gradient vector
is the canonical double commutator of the constrained Bridge Hamiltonian, or
must be modified by the reduced mass metric.
