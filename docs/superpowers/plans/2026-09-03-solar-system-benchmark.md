# Solar System Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reproducible 100-year full Solar System BridgeLF(4,2)/WHFast/IAS15 benchmark with Python plots and speedup/error summaries.

**Architecture:** A standalone C validation executable owns one deterministic fixture and runs each integrator, emitting sampled states/timings as CSV. A Python analysis script reads the CSV, computes IAS15-relative errors and speedups, and writes figures plus JSON.

**Tech Stack:** C99, REBOUND/rebound_bridge APIs, CMake validation target, Python 3 with numpy/matplotlib.

---

### Task 1: Add benchmark executable
**Files:** Create `validation/solar_system_full/solar_system_full_benchmark.c`; modify `CMakeLists.txt`.
- [ ] Define named bodies and deterministic masses/orbits, including dwarf planets and satellite groups.
- [ ] Implement three runs with identical initial states: IAS15, WHFast (`dt=0.25 d`), Bridge EOS LF(4,2) (`dt_outer=7 d`).
- [ ] Sample every configurable interval and write `method,time,body,x,y,z` rows plus elapsed seconds.
- [ ] Add CLI options `--years`, `--sample-years`, `--csv`, `--dt-whfast-days`, `--dt-bridge-days`.
- [ ] Register `solar_system_full_benchmark` under `REBOUND_BRIDGE_BUILD_VALIDATION`.

### Task 2: Add Python analysis and plotting
**Files:** Create `validation/solar_system_full/analyze_benchmark.py`.
- [ ] Load CSV, align method/body/time rows, and calculate relative Euclidean position error against IAS15.
- [ ] Emit `errors_timeseries.png`, `errors_summary.png`, `speedup.png`, and `summary.json`.
- [ ] Support `--csv`, `--outdir`, and `--title`; use a non-interactive matplotlib backend.

### Task 3: Document and verify
**Files:** Create `validation/solar_system_full/README.md`.
- [ ] Document build, 100-year run command, output schema, and non-fair timestep caveat.
- [ ] Build the target, run a short smoke test, then run the requested 100-year benchmark and Python analysis.
- [ ] Verify CSV row counts, finite errors, and positive timings; report generated artifact paths.
