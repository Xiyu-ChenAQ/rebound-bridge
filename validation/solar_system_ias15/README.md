# Solar System IAS15 Validation

This directory contains Solar System validation experiments comparing hierarchical Bridge integrations against IAS15.

The main C programs are:

- `compare_solar_system_ias15.c` - Bridge vs IAS15 comparison with orbital diagnostics and conservation metrics.
- `bench_whfast_bridge.c` - runtime comparison between fair global WHFast and Bridge.

Build and run an example case from the repository root:

```sh
cmake -S . -B build -DREBOUND_BRIDGE_BUILD_VALIDATION=ON
cmake --build build --config Release --target solar_system_ias15_compare
build/Release/solar_system_ias15_compare --years 100 --samples 101 --csv validation/solar_system_ias15/out/solar_system_ias15_compare.csv
```

On single-configuration generators such as Ninja or Unix Makefiles, the
executable may be at `build/solar_system_ias15_compare`.

Large generated output directories are not checked in.  For a compact first
review, use `../earth_moon_whfast_bridge/` for the WHFast/Bridge/IAS15 triangle
comparison.
