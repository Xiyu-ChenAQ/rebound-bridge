from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import matplotlib.pyplot as plt
import pytest


PAPER = Path(__file__).resolve().parents[1]


def load_scaling():
    spec = importlib.util.spec_from_file_location(
        "paper_runtime_scaling", PAPER / "scaling" / "plot_runtime_scaling.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_source_archive_has_no_required_git_metadata(tmp_path):
    completed = subprocess.run([
        sys.executable, "-c",
        "import json, runpy, sys; from pathlib import Path; "
        "script = runpy.run_path(sys.argv[1]); "
        "print(json.dumps(script['git_metadata'](Path(sys.argv[2]))))",
        str(PAPER / "multiple_subsystems" / "run_experiment.py"), str(tmp_path),
    ], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert json.loads(completed.stdout) == {"git_head": None, "git_status": None}


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def scaling_rows():
    return [
        {"n_outer": count, "bridge_method": "eos-lf4-2", "years": 100,
         "dt_outer_days": 4, "binary_ratio": 64, "whfast_dt_days": 0.0625,
         "repeats": 3, "warmups": 1, "whfast_mean_seconds": ratio,
         "bridge_mean_seconds": 1, "mean_speedup_whfast_over_bridge": ratio}
        for count, ratio in [(0, 2.5), (128, 30.0)]
    ]


def test_scaling_annotations_follow_new_measurements(tmp_path, monkeypatch):
    script = load_scaling()
    source = tmp_path / "new timings.csv"
    write_csv(source, scaling_rows())
    output = tmp_path / "figures" / "runtime_scaling_column"
    monkeypatch.setattr(sys, "argv", ["plot", "--csv", str(source), "--output", str(output)])
    # Inspect the real figure after writing its artifacts; only defer cleanup.
    close = script.plt.close
    monkeypatch.setattr(script.plt, "close", lambda *args: None)
    try:
        script.main()
        axis = script.plt.gcf().axes[0]
        labels = [text.get_text() for text in axis.texts]
        assert r"$2.50\times$" in labels
        assert r"$30.0\times$" in labels
        assert axis.get_ylim()[1] > 30
        assert output.with_suffix(".pdf").is_file()
        assert output.with_suffix(".png").is_file()
    finally:
        close("all")


def test_galilean_floor_annotations_follow_finest_kdk_sample(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "paper_convergence_figures", PAPER / "plot_convergence_figures.py"
    )
    script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(script)
    rows = [{"method": method, "sweep": sweep,
             "dt_outer_days": step if sweep == "outer" else 0.5,
             "dt_inner_days": step / 200 if sweep == "outer" else step / 10,
             "moon_error_km": step**order, "laplace_error_deg": 0.1 * step**order}
            for method, order in [("KDK / order2", 2), ("EOS LF4_2", 4)]
            for sweep in ["outer", "inner"] for step in [0.125, 0.25, 0.5]]
    rows.append({"method": "KDK / order2", "sweep": "inner",
                 "dt_outer_days": 0.5, "dt_inner_days": 0.002,
                 "moon_error_km": 37.0, "laplace_error_deg": 0.03})
    close = script.plt.close
    monkeypatch.setattr(script.plt, "close", lambda *args: None)
    try:
        script.galilean(list(reversed(rows)), tmp_path)
        for axis, value in zip(script.plt.gcf().axes[2:], [37.0, 0.03]):
            annotation = next(text for text in axis.texts
                              if text.get_text() == "Outer-coupling floor")
            assert annotation.xy == pytest.approx((0.002, value))
        assert (tmp_path / "galilean_two_timestep_convergence.pdf").is_file()
        assert (tmp_path / "galilean_two_timestep_convergence.png").is_file()
    finally:
        close("all")


@pytest.mark.parametrize("field,value", [
    ("bridge_mean_seconds", "nan"), ("whfast_mean_seconds", 0),
    ("binary_ratio", 32), ("bridge_method", "order2"),
])
def test_scaling_rejects_invalid_measurements(tmp_path, monkeypatch, field, value):
    script = load_scaling()
    rows = scaling_rows()
    rows[0][field] = value
    if field == "whfast_mean_seconds":
        rows[0]["mean_speedup_whfast_over_bridge"] = 0
    source = tmp_path / "bad.csv"
    write_csv(source, rows)
    output = tmp_path / "bad figure"
    monkeypatch.setattr(sys, "argv", ["plot", "--csv", str(source), "--output", str(output)])
    with pytest.raises(ValueError):
        script.main()
    assert not output.with_suffix(".pdf").exists()
    plt.close("all")


@pytest.mark.parametrize("selection,name,audit_key", [
    ("galilean", "galilean_two_timestep_convergence", "galilean_slopes"),
    ("generalized_order", "generalized_order_crossover_column", "crossover"),
    ("multiple_subsystems", "multiple_subsystems", "multiple"),
])
def test_selected_plot_uses_only_explicit_inputs(tmp_path, selection, name, audit_key):
    source = tmp_path / "input data"
    source.mkdir()
    output = tmp_path / "output figures"
    command = [sys.executable, str(PAPER / "plot_convergence_figures.py"),
               "--only", selection, "--output-dir", str(output),
               "--galilean-summary", str(source / "galilean.csv"),
               "--crossover-csv", str(source / "crossover.csv"),
               "--multiple-data", str(source / "multiple")]
    if selection == "galilean":
        rows = [{"method": method, "sweep": sweep,
                 "dt_outer_days": step if sweep == "outer" else 0.5,
                 "dt_inner_days": step / 200 if sweep == "outer" else step / 10,
                 "moon_error_km": step**order, "laplace_error_deg": 0.1 * step**order}
                for method, order in [("KDK / order2", 2), ("EOS LF4_2", 4)]
                for sweep in ["outer", "inner"] for step in [0.125, 0.25, 0.5]]
        write_csv(source / "galilean.csv", rows)
    elif selection == "generalized_order":
        rows = [{"method": "lf42", "coupling_scale": scale, "dt_days": step,
                 "max_position_error_au": 1e-6 * scale * step**2,
                 "reference_floor_au": 1e-15, "effective_epsilon": 1e-3 * scale}
                for scale in [0.01, 0.1, 1, 5] for step in [0.25, 0.5, 1]]
        write_csv(source / "crossover.csv", rows)
    else:
        multiple = source / "multiple"
        multiple.mkdir()
        rows = [{"coupling": coupling, "method": method, "epsilon": 1e-12,
                 "h_days": step, "error_full_model": step**2 * 1e-3,
                 "error_same_model": step**2 * 1e-4}
                for coupling in ["resolved", "monopole"] for method in ["lf4_2", "kdk"]
                for step in [0.03125, 0.25, 2]]
        (multiple / "summary.json").write_text(json.dumps({
            "metrics": rows, "default_vs_full_reference_model_gap": 1e-6,
            "reference_tolerance_check": {"resolved": 1e-14, "monopole": 1e-14},
            "delta_U_log_slope_last5": -5,
        }), encoding="utf-8")
        write_csv(multiple / "spatial_scaling.csv", [
            {"D_au": separation, "delta_U": 1e-8 / separation**5}
            for separation in [0.4, 0.8, 1.6, 3.2, 6.4]
        ])
    completed = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    audit_name = "plot_audit_" + selection + ".json"
    assert {path.name for path in output.iterdir()} == {
        name + ".pdf", name + ".png", audit_name
    }
    audit = json.loads((output / audit_name).read_text(encoding="utf-8"))
    assert set(audit) == {audit_key}
    if selection == "galilean":
        assert audit[audit_key]["KDK / order2:moon_error_km"] == pytest.approx(2)
        assert audit[audit_key]["EOS LF4_2:moon_error_km"] == pytest.approx(4)
