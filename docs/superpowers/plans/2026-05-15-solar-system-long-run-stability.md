# Solar System Long-Run Stability Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the existing standalone solar-system validation suite to run a 2000-year bridge-vs-IAS15 comparison and generate plots that expose long-term stability, bounded error, and residual growth.

**Architecture:** Keep the validation work isolated under `validation/solar_system_ias15/`. Reuse the current comparison executable as the single data producer, but make the duration and sampling configurable so the same tool can serve short comparisons and long stability runs. The Python plotting script should become contract-driven: it reads one CSV schema, then emits both summary plots and long-run stability plots from the same data.

**Tech Stack:** C99, REBOUND/REBOUND-bridge C API, Python 3, matplotlib, CMake, PowerShell.

---

### Task 1: Lock in the long-run plot contract

**Files:**
- Create: `validation/solar_system_ias15/test_plot_solar_system_ias15.py`
- Modify: `validation/solar_system_ias15/plot_solar_system_ias15.py`

- [ ] **Step 1: Write the failing test**

```python
def test_long_run_plot_contract(tmp_path):
    ...
    assert (tmp_path / "bridge_longterm_stability.png").exists()
    assert (tmp_path / "bridge_longterm_envelope.png").exists()
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m unittest validation.solar_system_ias15.test_plot_solar_system_ias15 -v`
Expected: FAIL because the new long-run plotting entry point does not exist yet.

- [ ] **Step 3: Write minimal plotting helpers**

Add `plot_long_term_stability(data, outdir)` and keep the existing comparison plots intact.

- [ ] **Step 4: Run the test to verify it passes**

Run: `python -m unittest validation.solar_system_ias15.test_plot_solar_system_ias15 -v`
Expected: PASS and the two long-run PNGs are created in the temp directory.

- [ ] **Step 5: Commit**

```powershell
git add validation/solar_system_ias15/test_plot_solar_system_ias15.py validation/solar_system_ias15/plot_solar_system_ias15.py
git commit -m "test: add long-run stability plot contract"
```

### Task 2: Parameterize the comparison executable for a 2000-year run

**Files:**
- Modify: `validation/solar_system_ias15/compare_solar_system_ias15.c`

- [ ] **Step 1: Write the failing build/run expectation**

Add CLI arguments or compile-time constants so the executable can be invoked as a 2000-year run without editing source each time. The new default run should target:

```text
years = 2000.0
samples = 2000
```

- [ ] **Step 2: Run the existing validation build**

Run: `cmake --build build-validation --config Release --target solar_system_ias15_compare`
Expected: build still passes after the new arguments are wired in.

- [ ] **Step 3: Implement the minimal changes**

Keep the existing CSV columns, but make the runtime duration configurable and preserve the current initial-condition setup and reference trajectory logic.

- [ ] **Step 4: Run the executable in long-run mode**

Run: `build-validation\Release\solar_system_ias15_compare.exe --years 2000 --samples 2000`
Expected: writes `validation/solar_system_ias15/out/solar_system_ias15_compare.csv` with 2000-year coverage.

- [ ] **Step 5: Commit**

```powershell
git add validation/solar_system_ias15/compare_solar_system_ias15.c
git commit -m "feat: add 2000-year validation mode"
```

### Task 3: Wire the long-run workflow and update docs

**Files:**
- Modify: `validation/solar_system_ias15/run.ps1`
- Modify: `README.md`

- [ ] **Step 1: Update the run script**

Make the script run the 2000-year mode and call the plotting script for both the existing comparison plots and the long-run stability plots.

- [ ] **Step 2: Update the README validation section**

Document that the validation suite now has a 2000-year long-run mode and list the generated figures in `validation/solar_system_ias15/out/`.

- [ ] **Step 3: Run the full workflow**

Run: `validation\solar_system_ias15\run.ps1`
Expected: CSV plus `bridge_vs_ias15_metrics.png`, `bridge_vs_ias15_residuals.png`, `bridge_host_errors.png`, `bridge_longterm_stability.png`, and `bridge_longterm_envelope.png`.

- [ ] **Step 4: Commit**

```powershell
git add validation/solar_system_ias15/run.ps1 README.md
git commit -m "docs: describe long-run solar-system validation"
```

### Task 4: Verify the long-run evidence

**Files:**
- Inspect: `validation/solar_system_ias15/out/solar_system_ias15_compare.csv`
- Inspect: generated PNGs under `validation/solar_system_ias15/out/`

- [ ] **Step 1: Re-run the validation workflow from scratch**

Run the PowerShell workflow again after a clean rebuild.

- [ ] **Step 2: Check the plots**

Confirm the long-run energy and angular-momentum curves stay bounded and that the envelope plots show no runaway drift over the 2000-year window.

- [ ] **Step 3: Report the result**

Summarize whether the bridge remains stable, where it deviates from IAS15, and whether the long-run evidence is sufficient for the repository README.

