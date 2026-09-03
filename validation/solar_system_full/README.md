# Full Solar System Benchmark

Run from the repository root:

```sh
python validation/solar_system_full/solar_system_full_benchmark.py --years 100 --sample-years 1 --csv validation/solar_system_full/results.csv
python validation/solar_system_full/analyze_benchmark.py --csv validation/solar_system_full/results.csv --outdir validation/solar_system_full/results
```

The production C benchmark uses the real EOS LF(4,2) implementation:

```sh
build-bench/Release/solar_system_whfast_bridge_bench.exe --years 100 --dt-outer-days 6.999808355691836 --earth-moon-ratio 20 --jovian-ratio 80 --bridge-method eos-lf4-2 --repeats 1 --warmups 0
```

This measured Bridge `0.353042600 s` versus WHFast `0.975974600 s` (`2.764466951x`). The Python file is only a fixture/plot prototype and must not be used to claim EOS LF(4,2) accuracy; a state-output C driver is still needed for trajectory-error plots.
