# Earth-Moon WHFast/Bridge validation

This validation compares the hierarchical Bridge split against a naive global WHFast baseline and IAS15 for a real multi-satellite Solar System setup.

The system contains:

- Sun and major planets.
- Earth-Moon as a resolved Bridge subsystem.
- Jupiter-Io-Europa-Ganymede as a resolved Bridge subsystem.

The purpose is to test whether a hierarchical Bridge factorization is better matched to multi-center satellite systems than a single global WHFast integration.

## Reproduce

From the repository root:

```sh
cmake -S . -B build -DREBOUND_BRIDGE_BUILD_VALIDATION=ON
cmake --build build --config Release --target earth_moon_whfast_bridge_compare
```

Then run the executable directly. Example:

```sh
build/Release/earth_moon_whfast_bridge_compare --years 100 --samples 101 --dt-outer-days 0.25 --earth-moon-ratio 20 --jovian-ratio 80 --bridge-method eos-lf4-2
```

On single-configuration generators such as Ninja or Unix Makefiles, the
executable may be at `build/earth_moon_whfast_bridge_compare`.

The checked-in `results/` directory contains the current CSV and PNG summary.

## Current 100-year result

The checked-in `results/` directory contains the latest local measurements used for the current paper draft.

For `dt_outer = 0.25 d`, `earth_moon_ratio = 20`, `jovian_ratio = 80`:

- Bridge vs IAS15 Moon relative error: about `28.8 km`.
- Global WHFast vs IAS15 Moon relative error: about `691.5 km`.
- Bridge vs global WHFast Moon relative difference: about `720.3 km`.

This does not mean Bridge is generally more accurate than WHFast. It means the hierarchical Bridge split is better matched to this multi-center satellite problem than a naive global WHFast baseline.
