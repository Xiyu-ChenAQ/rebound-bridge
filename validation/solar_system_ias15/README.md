# Solar System IAS15 Validation

This directory contains Solar System validation experiments comparing hierarchical Bridge integrations against IAS15.

The main C programs are:

- `compare_solar_system_ias15.c` - Bridge vs IAS15 comparison with orbital diagnostics and conservation metrics.
- `bench_whfast_bridge.c` - runtime comparison between fair global WHFast and Bridge.

Build and run an example case from the repository root:

```sh
cmake -S . -B build
cmake --build build --config Release --target solar_system_ias15_compare
build/Release/solar_system_ias15_compare --years 100 --samples 101 --csv validation/solar_system_ias15/out/solar_system_ias15_compare.csv
```

Generated output directories such as `out_2000yr/`, `out_20000yr/`, and `benchmark_figures/` are historical or regenerated results. They are not the recommended first entry point for reviewers; use `../earth_moon_whfast_bridge/` for the compact WHFast/Bridge/IAS15 triangle comparison.
