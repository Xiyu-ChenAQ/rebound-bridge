# Pluto-Charon Binary Validation

This directory validates hierarchical Bridge behavior on a fast binary subsystem with a large mass ratio compared with ordinary planet-satellite systems.

The main C programs are:

- `compare_pluto_charon.c` - Bridge vs IAS15 comparison for Pluto-Charon-like binary cases.
- `bench_pluto_charon.c` - fair WHFast/Bridge runtime benchmark.

Build and run an example case from the repository root:

```sh
cmake -S . -B build
cmake --build build --config Release --target pluto_charon_compare
build/Release/pluto_charon_compare --case outer --years 100 --samples 101 --dt-outer-days 2 --binary-ratio 64
```

Plotting:

- `plot_pluto_charon.py` - plots accuracy, phase, COM, and runtime diagnostics from sweep CSV files.

This suite is useful for testing whether Bridge remains stable for hierarchical binary systems where the subsystem barycenter is not close to the primary body.
