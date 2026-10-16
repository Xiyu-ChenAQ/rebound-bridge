# Reproducing the paper

These instructions cover the `1.0.0b1` beta release.

This directory contains the experiment programs, initial conditions, analysis
and plotting scripts for the seven figures in the REBOUND_BRIDGE manuscript.
It builds against the library and REBOUND sources included in the repository
or source distribution. It does not require the manuscript or the author's
local `research/` directory.

The Python wheel installs the usable library. For these experiments, use the
**source archive** (`rebound_bridge-1.0.0b1.tar.gz`) or a complete repository
checkout. Installing the wheel alone does not install `tests/paper/`.

## Start here

Use Python 3.10 or newer, CMake 3.20 or newer, and a C99 compiler. On Windows,
run the build in a Visual Studio developer terminal with CMake available.
From the extracted source root:

```sh
python -m pip install -r tests/paper/requirements.txt
python tests/paper/run_all.py --unit-tests-only
python tests/paper/run_all.py --mode smoke --output-dir build/paper-smoke
```

The smoke run builds all nine C programs, runs the registered CTest checks, and
exercises all experiment families. Most trajectories last eight days; the
generalized-order experiment retains its inexpensive 128-day configuration,
and the two-subsystem experiment uses its reduced step grid. These checks
verify execution and output, not the long-duration results reported in the
paper. Smoke outputs are not named or plotted as full manuscript results.

For the full set of manuscript experiments:

```sh
python tests/paper/run_all.py --mode full --dry-run
python tests/paper/run_all.py --mode full --output-dir build/paper-full
```

`--dry-run` prints a JSON command plan and creates no files. `full` includes
the 20-year repeated speed tests and four 10,000-year realizations of each
long-term method, so it can take a long time. Run performance comparisons on
an otherwise idle machine. Runtime ratios depend on the hardware and compiler.

The result directory must be fresh. A normal run refuses to overwrite an
existing nonempty directory. Use a separate output directory to repeat a run.
Every run records commands, return codes, logs, source and executable hashes,
output hashes, platform, Python version, and completion status in
`run_manifest.json`. A failed command stops the run and leaves its log and
partial output available for inspection.

## Figures and experiments

| Figure / section | `--experiments` value | Experiment and settings | Final PDF in `figures/` |
| --- | --- | --- | --- |
| Fig. 1 / 5.1 | `convergence` | Galilean system, 4 yr; outer steps 0.125, 0.25, 0.5, 1 d with ratio 200; inner ratios 5–320 at outer step 0.5 d | `galilean_two_timestep_convergence.pdf` |
| Fig. 2 / 5.2 | `generalized_order` | Synthetic star–planet–moon system, 128 d; residual scale 0.01–5, outer steps 0.25–8 d; IAS15 components; fine references at 1/32 and 1/64 d | `generalized_order_crossover_column.pdf` |
| Fig. 3 / 5.3 | `speed_accuracy` | Galilean system, 20 yr; global SABA and WHFast baselines; KDK/LF(4,2) at 0.125–1 d, inner step 0.0125 d, five interleaved repetitions with 5-d output | `galilean_bridge_vs_saba_whfast_20yr_exactcadence_column.pdf` |
| Fig. 4 / 5.4 | `scaling` | Sun–binary system with 0, 4, 8, 16, 32, 64, 128 added bodies; nominal 100 yr; outer step 4 d, inner/global step 0.0625 d; three runs after one warm-up | `runtime_scaling_column.pdf` |
| Fig. 5 / 5.5 | `multiple_subsystems` | Two resolved subsystems, 128 d; default/full coupling, independent IAS15 references, outer steps 2 to 1/256 d; fixed-geometry potential evaluated with 70-digit arithmetic | `multiple_subsystems.pdf` |
| Fig. 6 / 5.6 | `long_energy` | Four perturbed Galilean realizations, 10,000 yr; Bridge outer/inner steps 0.5/0.0125 d; corrected WHFast step 0.0085 d | `galilean_long_ensemble_10000yr_r4_column.pdf` |
| Fig. 7 / 6 | `exomoon` | Planet at 0.55 au; 20 yr at three Bridge steps and three initial phases, direct IAS15, plus 100-yr Bridge/IAS15 comparison | `laplace_circulation_validation_column.pdf` |
| Section 4 | `reversibility` | Galilean system, 1,000 steps each way; KDK and LF(4,2), outer/inner steps 0.5/0.0125 d | `reversibility.csv` (no figure) |

In Fig. 3, all methods use the same 5-day output intervals over 20 years
(1,461 outputs), and each reported maximum is taken over this common grid.
The global SABA/WHFast comparisons remain single runs; KDK/LF(4,2) use five
interleaved repetitions. This `1.0.0b1` configuration replaces the earlier global
baseline schedule of ten outputs at 2, 4, ..., 20 yr. Its errors and runtimes
must therefore be recomputed rather than combined with the older results.

For example, run only the two-subsystem experiment:

```sh
python tests/paper/run_all.py --mode full --experiments multiple_subsystems --output-dir build/paper-two-subsystems
```

Several experiment names may follow `--experiments`. `--build-dir` selects
the CMake output directory (default `build/paper`); `--config` defaults to
`Release`. Use `--skip-build` only when those executables already match the
source version being tested.

## Redrawing existing results

```sh
python tests/paper/run_all.py --mode full --plots-only --output-dir build/paper-full
python tests/paper/run_all.py --mode full --plots-only --experiments generalized_order --output-dir build/paper-full
```

This only runs plotting commands; it does not rebuild or integrate. It writes
a separate timestamped plot manifest and leaves `run_manifest.json` intact.
Figure 1 can be redrawn from its summary CSV without retaining all 22 raw
trajectories. The current layouts for Figs. 1, 2 and 5 are generated by
`plot_convergence_figures.py`; the two-panel image produced internally by the
two-subsystem runner is an additional diagnostic, not the manuscript figure.

## Numerical settings and version differences

The fixed Galilean states in `common/galilean_horizons_state.h` are from
NASA/JPL Horizons at 2026 July 8, 00:00 TDB. No network query is needed.
The other synthetic initial conditions and deterministic perturbation seeds
are included in the experiment sources. Units are au, yr and solar masses.

The library requires complete outer steps. Smoke output times are aligned
accordingly. The Horizons and exomoon drivers validate the integer output grid
before integration, advance a fixed number of outer steps per output, and
record the actual integration clock. In the Horizons comparison, the direct
reference is evaluated at that same clock. Non-divisible durations or output
counts are rejected rather than rounded to another grid.

For the nominal 100-year scaling test, 36,525 days is not a
multiple of the 4-day step. Both methods therefore stop after 9,131 complete
outer steps (36,524 days); the CSV records requested and actual durations
and step counts. Long-term integrations likewise use whole fixed steps, as
recorded in their trajectory times.

All Bridge benchmarks use the `1.0.0b1` default WHFast Jacobi-state reuse and direct
Jacobi velocity updates when the component configuration supports them.
Earlier convergence/speed, scaling and exomoon runs used Cartesian increments.
This changes coordinate-update arithmetic, not the Hamiltonian split or
composition coefficients, but new results must not be mixed with old runs.

The exomoon CSV retains `inner_resonant_angle_deg` (the angle using Europa's
longitude of periapse) and adds `theta1_deg` (the angle using Io's longitude
of periapse). They are different diagnostics. A missing column in an older
CSV must not be interpreted as evidence about the other angle.

Scientific reproduction should compare diagnostics and convergence under the
specified settings. Floating-point rounding and measured runtimes need not
match the manuscript bit for bit. Short tests, successful builds and redraws
alone do not establish a complete rerun of the long-duration experiments.

## Two-subsystem experiment

Figure 5 uses a synthetic five-body system to distinguish integration error
from the spatial approximation in default coupling. It tests both KDK and
LF(4,2) with two resolved subsystems; it is not a model of an observed system
or a runtime benchmark. The checks are:

1. Compare the residual kick with the derivative of each coupling potential,
   using independent analytic C forces and 70-digit Python differentiation.
2. Check that fully resolved coupling approaches direct Newtonian IAS15 as
   the outer step decreases.
3. Check that default coupling approaches its own Hamiltonian,
   `U_s0 + U_0t - U_00`. Its independent IAS15 reference distributes proxy
   reactions through the centre-of-mass Jacobian.
4. Hold the internal geometry fixed and vary the subsystem separation to
   measure the omitted potential's predicted `D^-5` scaling.

### Initial conditions and diagnostics

Units are au, yr and solar masses, with `G = 4*pi^2` and a unit-mass star.
The moon masses make the inter-subsystem approximation visible in a short run.

| Parameter | Subsystem 1 | Subsystem 2 |
| --- | ---: | ---: |
| Host mass | 0.001 | 0.0008 |
| Moon mass | 0.00005 | 0.00004 |
| Initial COM radius (au) | 1 | 1.45 |
| Internal semimajor axis (au) | 0.02 | 0.025 |
| Internal eccentricity | 0.05 | 0.04 |
| Internal inclination (rad) | 0.1 | -0.12 |
| Internal true anomaly (rad) | 0.3 | 1.8 |
| Internal node, periapsis argument (rad) | 0, 0 | 0, 0 |

Both subsystem COMs start on the positive x axis, with tangential speeds
`sqrt(G*(M_star + M_subsystem)/radius)`. The states are shifted to the full-system
barycentric frame and saved in every trajectory. Runs last 128 days, sampled
every 2 days including the initial state (65 epochs). The outer steps are
2, 1, 1/2, ..., 1/256 days, so every output lies on each tested step grid.

Components use adaptive IAS15 with tolerance `1e-12` to suppress component
error; their flows are not claimed to be exactly symplectic. The nominal inner
interval is `h/8`, within which IAS15 chooses its steps. Independent references
use tolerances `1e-13` and `1e-14`. LF(4,2) is also run with component tolerance
`1e-14` at `h = 1/8` and `1/256` days in both coupling modes.

The trajectory error is the maximum over all epochs and both subsystems of
the host-relative Cartesian position difference, divided by that subsystem's
initial semimajor axis. Momentum and angular momentum use reconstructed
physical particles; energy uses the Hamiltonian for the selected coupling
mode. These are finite-duration diagnostics, not long-term bounds.

The separate potential sweep uses `D = 0.4*1.5^k` au for `k = 0,...,9` at fixed
internal masses, offsets and orientations. The fit uses the largest five
separations. Subtracting the potentials at 70 decimal digits avoids cancellation
of the small omitted term; this is a spatial scaling check, not a trajectory fit.

### Running the experiment directly

The [experiment selector](#figures-and-experiments) runs this test as part of
the common reproduction workflow. To use its standalone runner instead:

```sh
cmake -S tests/paper -B build/paper
cmake --build build/paper --config Release --target multiple_subsystems
python tests/paper/multiple_subsystems/run_experiment.py --exe build/paper/Release/multiple_subsystems.exe --output build/paper-two-subsystems-standalone
```

For a single-configuration Unix build, use `build/paper/multiple_subsystems`
as the executable. Runs are sequential, with a 120-second timeout per child
execution, and require a fresh output directory. Add `--pilot` to use three
outer steps and skip component-tightening runs. To reanalyse saved trajectories,
use `--analyse-only --output <result-directory>` without `--exe`.

Outputs include raw trajectories, `kick_check.csv`, `convergence.csv`,
`local_orders.csv`, `spatial_scaling.csv`, `summary.json`, and a diagnostic
PDF/PNG. All four method/coupling sweeps, including default KDK, are retained.
The adjacent-step orders are recorded without assuming LF(4,2) remains fourth
order: its term proportional to `epsilon^2*h^2` can dominate. For the manuscript
figure layout, use the [redrawing workflow](#redrawing-existing-results).

## Source-package verification

From a checkout with the Python build dependencies installed:

```sh
python -m build --sdist
python tests/paper/checks/check_paper_sdist.py dist/rebound_bridge-1.0.0b1.tar.gz
```

The archive checker verifies that the current experiment and CMake source
files are present and match the checkout, and rejects generated results,
caches and binaries. Raw outputs are generated locally and are not bundled
in the source archive. Preserve the archive's SHA-256 alongside the run
manifest when sharing a reproduction.
