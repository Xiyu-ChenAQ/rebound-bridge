# Full Solar System Benchmark

Run from the repository root:

```sh
python validation/solar_system_full/solar_system_full_benchmark.py --years 100 --sample-years 1 --csv validation/solar_system_full/results.csv
python validation/solar_system_full/analyze_benchmark.py --csv validation/solar_system_full/results.csv --outdir validation/solar_system_full/results
```

WARNING: `solar_system_full_benchmark.py` is not a valid EOS LF(4,2) result. It is retained only as a fixture/plot prototype and uses REBOUND's leapfrog integrator as a placeholder. Do not cite its Bridge numbers. The production benchmark must call the C API with `REB_BRIDGE_METHOD_EOS_LF4_2`; that implementation is not yet wired to this extended fixture.
