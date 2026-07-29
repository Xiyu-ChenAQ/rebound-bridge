# Paper validation suite

This directory consolidates the numerical programs and plotting scripts used
for the five figures and the reversibility check in the MNRAS manuscript. It
is self-contained with respect to experiment code, but links against the
parent `rebound_bridge` source tree. Existing files under `MNRAS/` are not used
at runtime and are left unchanged.

## Contents

- `galilean/`: timestep convergence, speed--accuracy, long-term conservation,
  and reversibility experiments.
- `generalized_order/`: generalized-order crossover experiment.
- `scaling/`: global-WHFast/Bridge runtime scaling with slow outer bodies.
- `exomoon/`: Laplace-circulation validation, including Bridge and direct
  IAS15 modes.
- `common/`: the JPL Horizons initial state shared by the Galilean programs.
- `tests/`: Python analysis tests and the three automatic limiting-case checks.
- `results/`: generated CSV tables and figures. Large historical results are
  intentionally not copied.
- `run_all.py`: one build and execution entry point.

## Setup and smoke test

From the repository root:

```powershell
python -m pip install -r paper_test/requirements.txt
python paper_test/run_all.py
```

The default smoke mode builds all eight C executables, runs short instances of
each experiment (including Bridge and IAS15 exomoon cases), and runs the Python
tests. Output is written below `paper_test/results/smoke/`.

Useful commands:

```powershell
python paper_test/run_all.py --unit-tests-only
python paper_test/run_all.py --mode smoke --skip-build
python paper_test/run_all.py --mode full --dry-run
python paper_test/run_all.py --mode full
python paper_test/run_all.py --mode full --plots-only
```

`--dry-run` prints the exact experiment manifest without running it.
`--plots-only` regenerates figures from existing full-mode CSV files.

## Full reproduction

Full mode reproduces the inputs used for the manuscript's five figures:

1. Four-year two-timestep convergence sweeps.
2. The generalized-order crossover.
3. Twenty-year SABA, WHFast, corrected-WHFast, KDK, and LF(4,2)
   speed--accuracy measurements; KDK/LF(4,2) timings use five interleaved runs.
4. Same-session `1000 yr` cost calibration and four deterministic `10000 yr`
   Bridge and corrected-WHFast realizations.
5. The `0.55 au` exomoon tests at three phases, three Bridge outer steps,
   direct IAS15, and the additional `100 yr` Bridge/IAS15 check.
6. The independent 0/32/64/128-outer-body runtime scaling quoted in the text.

Full mode is intentionally expensive. Timing results depend on the machine,
compiler, and host load. In particular, the five-repeat 20-year benchmark and
the four-realization `10000 yr` ensemble should be run on an otherwise idle
machine. Generated CSV files can be large and are not source artifacts.
