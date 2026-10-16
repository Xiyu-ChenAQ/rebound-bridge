from __future__ import annotations

import csv
import hashlib
from pathlib import Path
import subprocess
import sys

import matplotlib.image as mpimg
import pytest


PAPER = Path(__file__).resolve().parents[1]


def write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_script(script, *arguments, cwd):
    completed = subprocess.run(
        [sys.executable, str(PAPER / script), *map(str, arguments)],
        cwd=cwd, capture_output=True, text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


def hashes(directory):
    return {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in directory.glob("*.csv")}


def assert_two_layouts(output, stem):
    for name in (stem, stem + "_column"):
        for suffix in (".pdf", ".png"):
            assert (output / (name + suffix)).stat().st_size > 1000
    horizontal = mpimg.imread(output / (stem + ".png"))
    column = mpimg.imread(output / (stem + "_column.png"))
    assert horizontal.shape[1] > horizontal.shape[0]
    assert column.shape[0] > column.shape[1]


def test_two_timestep_skip_run_reads_all_cases_in_output_dir(tmp_path):
    output = tmp_path / "custom convergence"
    for method in ("order2", "eos_lf4_2"):
        cases = [("outer", step, 200) for step in (0.125, 0.25, 0.5, 1)]
        cases += [("inner", 0.5, ratio) for ratio in (5, 10, 20, 40, 80, 160, 320)]
        for sweep, step, ratio in cases:
            token = f"{step:g}".replace(".", "p")
            write_rows(output / f"galilean_{sweep}_{method}_4yr_dt{token}_ir{ratio}.csv", [{
                "max_moon_error_au": step**2 * 1e-6, "laplace_diff_deg": step**2 * 1e-3,
                "bridge_energy_rel": -1e-10, "cumulative_speedup": 2,
            }])
    original = hashes(output)
    run_script("galilean/plot_two_timestep_convergence.py", "--skip-run",
               "--output-dir", output, cwd=tmp_path)
    summary = output / "galilean_two_timestep_convergence_summary.csv"
    with summary.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 22
    assert float(rows[0]["moon_error_km"]) == pytest.approx(2.3374667296875)
    assert float(rows[0]["energy_rel_abs"]) == 1e-10
    assert float(rows[0]["dt_inner_days"]) == 0.000625
    assert all(hashes(output)[name] == value for name, value in original.items())
    assert (output / "galilean_two_timestep_convergence.pdf").is_file()
    assert (output / "galilean_two_timestep_convergence.png").is_file()


def test_speed_accuracy_cli_keeps_summary_and_both_layouts_in_output(tmp_path):
    data = tmp_path / "global runs"
    cost = tmp_path / "bridge runs"
    output = tmp_path / "speed figures"
    cases = [("saba", ("0p1", "0p05", "0p025", "0p02", "0p01"))]
    cases += [(method, ("0p02", "0p01", "0p005", "0p002", "0p001"))
              for method in ("whfast", "whfast_corr11c2")]
    for method, tokens in cases:
        for token in tokens:
            write_rows(data / f"galilean_{method}_20yr_dt{token}_exactcadence.csv", [
                {"time_yr": 10, "cumulative_speedup": 1.5, "max_moon_error_au": 2e-6,
                 "laplace_diff_deg": 2e-3, "split_escaped": 0},
                {"time_yr": 20, "cumulative_speedup": 2, "max_moon_error_au": 1e-6,
                 "laplace_diff_deg": 1e-3, "split_escaped": 0},
            ])
    write_rows(cost / "kdk_lf42_cost_accuracy_summary.csv", [
        {"method": method, "method_label": method, "dt_outer_days": step,
         "speedup_median": step * 4, "speedup_q25": step * 3.5,
         "speedup_q75": step * 4.5, "max_moon_error_au": step**2 * 1e-6,
         "max_laplace_diff_deg": step**2 * 1e-3}
        for method in ("order2", "eos_lf4_2") for step in (0.125, 0.25, 0.5, 1)
    ])
    original_data, original_cost = hashes(data), hashes(cost)
    run_script("galilean/plot_20yr_speed_accuracy.py", "--data-dir", data,
               "--cost-dir", cost, "--output-dir", output, cwd=tmp_path)
    assert hashes(data) == original_data
    assert hashes(cost) == original_cost
    with (output / "galilean_global_baselines_20yr_exactcadence_summary.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 15
    assert float(rows[0]["max_moon_error_km"]) == pytest.approx(299.1957414)
    assert float(rows[0]["final_moon_error_km"]) == pytest.approx(149.5978707)
    assert float(rows[0]["max_laplace_diff_deg"]) == 2e-3
    assert_two_layouts(output, "galilean_bridge_vs_saba_whfast_20yr_exactcadence")


def test_exomoon_cli_preserves_unwrapping_and_produces_column_layout(tmp_path):
    data = tmp_path / "exomoon runs"
    output = tmp_path / "exomoon figures"
    names = ["a0p55_phase0_20yr_ias15.csv"]
    names += [f"a0p55_phase0_20yr_bridge_dt{token}.csv" for token in ("0p5", "0p25", "0p125")]
    names += [f"a0p55_phase{phase}_20yr_bridge_dt0p25.csv" for phase in (120, 240)]
    for name in names:
        write_rows(data / name, [
            {"time_yr": time, "laplace_deg": angle % 360, "ever_escaped": 0}
            for time, angle in [(0, 0), (1, -120), (2, -240), (3, -360),
                                (4, -480), (5, -600), (20, -720)]
        ])
    original = hashes(data)
    run_script("exomoon/plot_laplace_circulation_validation.py", "--data-dir", data,
               "--output-dir", output, cwd=tmp_path)
    assert hashes(data) == original
    with (output / "laplace_circulation_validation_summary.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 6
    assert float(rows[0]["laplace_span_deg"]) == pytest.approx(720)
    assert float(rows[0]["laplace_net_turns"]) == pytest.approx(-2)
    assert float(rows[0]["first_full_circulation_yr"]) == 3
    assert_two_layouts(output, "laplace_circulation_validation")
