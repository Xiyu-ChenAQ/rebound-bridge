# Source Layout

The `rebound_bridge` source tree is split by bridge responsibility rather than by numerical method alone.

Recommended reading order:

1. `bridge.c` - public API, lifecycle, subsystem registration, method selection, and workspace ownership.
2. `bridge_advance.c` - outer advance driver and time synchronization.
3. `bridge_kick.c` - cross-kick operators coupling main and subsystem simulations.
4. `bridge_eos.c` - EOS LF composition kernels.
5. `bridge_internal.h` - private structures, inline utilities, and internal function declarations.

The public API is declared in `src/include/rebound_bridge.h`.
