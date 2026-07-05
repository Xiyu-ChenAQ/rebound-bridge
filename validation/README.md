# Validation

This directory contains standalone validation programs for `rebound_bridge`.

The validation programs compare hierarchical Bridge integrations against IAS15 reference integrations and fair global WHFast baselines. They are built as normal CMake targets; no wrapper scripts are required.

## Suites

- `earth_moon_whfast_bridge/` compares hierarchical Bridge, global WHFast, and IAS15 on a real Solar System setup with resolved Earth-Moon and Galilean satellite subsystems.
- `solar_system_ias15/` contains longer Solar System Bridge-vs-IAS15 validation programs.
- `pluto_charon_binary/` validates fast binary subsystems using a Pluto-Charon-like configuration.
- `validation_timer.h` provides the shared wall-clock timer used by benchmark programs.

## Build

From the repository root:

```sh
cmake -S . -B build
cmake --build build --config Release
```

Useful targets:

- `earth_moon_whfast_bridge_compare`
- `solar_system_ias15_compare`
- `solar_system_whfast_bridge_bench`
- `pluto_charon_compare`
- `pluto_charon_whfast_bridge_bench`

## Recommended Entry Point

Start with:

```sh
cmake --build build --config Release --target earth_moon_whfast_bridge_compare
build/Release/earth_moon_whfast_bridge_compare --years 100 --samples 101 --dt-outer-days 0.25 --earth-moon-ratio 20 --jovian-ratio 80 --bridge-method eos-lf4-2
```

On single-configuration generators, the executable may be at `build/earth_moon_whfast_bridge_compare`.

The checked-in summary figure is:

```text
earth_moon_whfast_bridge/results/earth_moon_whfast_bridge_triangle.png
```

It shows that, for this multi-center satellite problem, the hierarchical Bridge split is closer to IAS15 than a naive global WHFast baseline at the same minimum timestep.

## Output Policy

Small CSV/PNG summaries may be checked in when they document a specific comparison. Large long-term output directories should remain ignored by git.
