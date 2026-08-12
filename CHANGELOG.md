# Changelog

All notable changes to `rebound_bridge` are documented in this file.

## [0.1.0] - 2026-08-12

Initial research release.

### Added

- A C99 Bridge scheduler for a main REBOUND simulation and resolved fast
  subsystems with independent timesteps and integrators.
- Momentum-conserving resolved-minus-monopole residual kicks for fixed
  hierarchical systems.
- Stable second-order KDK and LF(4,2) Bridge compositions.
- Optional resolved coupling between registered subsystems.
- CMake package installation and exported targets for downstream projects.
- Example, validation, regression, reversibility, and paper-reproduction
  suites.

### Experimental

- LF(8,6,4), force-gradient, and commutator-corrector interfaces remain
  validation targets and are not recommended as production defaults.

[0.1.0]: https://github.com/Xiyu-ChenAQ/rebound-bridge/releases/tag/v0.1.0
