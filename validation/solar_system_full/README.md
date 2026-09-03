# Full Solar System Benchmark

Run from the repository root:

```sh
python validation/solar_system_full/solar_system_full_benchmark.py --years 100 --sample-years 1 --csv validation/solar_system_full/results.csv
python validation/solar_system_full/analyze_benchmark.py --csv validation/solar_system_full/results.csv --outdir validation/solar_system_full/results
```

The fixture includes the Sun, eight planets, Ceres, Pluto, Eris, Haumea, Makemake, the Moon, Phobos/Deimos, four Galilean moons, and four major Saturnian moons. IAS15 is the trajectory reference. BridgeLF(4,2) is represented with the requested 7-day outer cadence; WHFast uses a conservative 0.25-day global step, so the timing comparison is intentionally not cost-balanced.
