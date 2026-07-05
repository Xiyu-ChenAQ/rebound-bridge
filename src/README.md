# Source Layout

The `rebound_bridge` source tree is split by bridge responsibility rather than by numerical method alone.

Recommended reading order:

1. `bridge.c` - public API, lifecycle, subsystem registration, method selection, and workspace ownership.
2. `bridge_advance.c` - outer advance driver and time synchronization.
3. `bridge_kick.c` - cross-kick operators coupling main and subsystem simulations.
4. `bridge_eos.c` - EOS LF composition kernels.
5. `bridge_fg4_experimental.c` - experimental forward force-gradient composition.
6. `bridge_force_gradient_experimental.c` - force-gradient correction helpers and diagnostics.
7. `bridge_processor.c` - experimental commutator processor.
8. `bridge_examples.c` - reusable example-system constructors.
9. `bridge_internal.h` - private structures, inline utilities, and internal function declarations.

The stable public API is declared in `include/rebound_bridge.h`. Experimental options are declared in `include/rebound_bridge_experimental.h`.
