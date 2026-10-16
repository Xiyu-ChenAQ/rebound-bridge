from __future__ import annotations

from pathlib import Path
import json
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import run_all


def test_exomoon_paper_suite_contains_direct_and_long_checks() -> None:
    command = run_all.exomoon_command(Path("experiment"), "ias15", 100, 0.5, 0, Path("out.csv"))

    assert "ias15" in command
    assert command[command.index("--samples") + 1] == "73050"
    assert command[-1] == Path("out.csv")


def test_default_mode_is_fast_smoke() -> None:
    args = run_all.parse_args([])

    assert args.mode == "smoke"
    assert not args.plots_only


def test_all_paper_executables_are_in_build_manifest() -> None:
    assert "complex_speed_ratio" in run_all.TARGETS
    assert "bridge_limiting_cases_test" in run_all.TARGETS
    assert "galilean_exomoon_survival" in run_all.TARGETS
    assert "multiple_subsystems" in run_all.TARGETS


def dry_run(output: Path, *options: str) -> dict:
    completed = subprocess.run(
        [sys.executable, str(Path(run_all.__file__)), "--dry-run",
         "--output-dir", str(output), *options],
        capture_output=True, text=True, check=True,
    )
    return json.loads(completed.stdout)


def test_full_plan_covers_seven_figures_without_creating_results(tmp_path: Path) -> None:
    output = tmp_path / "new results"
    plan = dry_run(output, "--mode", "full")
    assert not output.exists()
    assert len(plan["figures"]) == 7
    commands = [row["argv"] for row in plan["commands"]]
    assert any("--n-outer" in command and "128" in command for command in commands)
    assert any(any("run_experiment.py" in arg for arg in command) for command in commands)
    assert any(any("plot_convergence_figures.py" in arg for arg in command) for command in commands)


def test_smoke_output_times_are_on_the_outer_grid(tmp_path: Path) -> None:
    plan = dry_run(tmp_path / "smoke", "--skip-build")
    for row in plan["commands"]:
        command = row["argv"]
        if "--years" not in command or "--dt-outer-days" not in command:
            continue
        days = float(command[command.index("--years") + 1]) * 365.25
        step = float(command[command.index("--dt-outer-days") + 1])
        samples = int(command[command.index("--samples") + 1]) if "--samples" in command else 1
        ratio = days / samples / step
        assert abs(ratio - round(ratio)) < 1e-10, command


def test_full_speed_comparison_uses_common_five_day_outputs(tmp_path: Path) -> None:
    plan = dry_run(tmp_path / "full", "--mode", "full", "--experiments", "speed_accuracy")
    for row in plan["commands"]:
        command = row["argv"]
        if row["experiment"] != "speed_accuracy" or row["kind"] != "integration":
            continue
        days = float(command[command.index("--years") + 1]) * 365.25
        if "--samples" in command:
            samples = int(command[command.index("--samples") + 1])
            assert days / samples == 5.0
            step = float(command[command.index("--dt-days") + 1])
            assert abs(5.0 / step - round(5.0 / step)) < 1e-10
        else:
            assert command[command.index("--output-interval-days") + 1] == "5"


def test_selected_plot_plan_does_not_schedule_integrations(tmp_path: Path) -> None:
    plan = dry_run(tmp_path / "data", "--mode", "full", "--plots-only",
                   "--experiments", "generalized_order")
    assert len(plan["figures"]) == 1
    assert all(row["kind"] == "plot" for row in plan["commands"])
    assert "generalized_order_crossover_column.pdf" in plan["figures"][0]


def test_existing_results_are_not_overwritten(tmp_path: Path) -> None:
    output = tmp_path / "existing"
    output.mkdir()
    sentinel = output / "valuable.csv"
    sentinel.write_text("original\n", encoding="utf-8")
    completed = subprocess.run(
        [sys.executable, str(Path(run_all.__file__)), "--skip-build",
         "--output-dir", str(output), "--experiments", "reversibility"],
        capture_output=True, text=True,
    )
    assert completed.returncode != 0
    assert "not empty" in completed.stderr
    assert sentinel.read_text(encoding="utf-8") == "original\n"


def test_smoke_run_cannot_be_redrawn_as_manuscript_results(tmp_path: Path) -> None:
    output = tmp_path / "smoke"
    output.mkdir()
    (output / "run_manifest.json").write_text('{"mode": "smoke"}', encoding="utf-8")
    completed = subprocess.run(
        [sys.executable, str(Path(run_all.__file__)), "--mode", "full", "--plots-only",
         "--output-dir", str(output), "--experiments", "generalized_order"],
        capture_output=True, text=True,
    )
    assert completed.returncode != 0
    assert "smoke results" in completed.stderr
    assert not (output / "figures").exists()


def test_commands_use_actual_build_layout_after_configuration(tmp_path: Path) -> None:
    suffix = ".exe" if sys.platform == "win32" else ""
    name = "galilean_horizons_compare" + suffix
    for actual, guessed in ((tmp_path / name, tmp_path / "Release" / name),
                             (tmp_path / "Release" / name, tmp_path / name)):
        actual.parent.mkdir(parents=True, exist_ok=True)
        actual.write_bytes(b"executable fixture")
        command = [sys.executable, "analysis.py", "--exe", str(guessed), "--csv", "out.csv"]
        resolved = run_all.resolve_executables(command, tmp_path, "Release")
        assert resolved == [sys.executable, "analysis.py", "--exe", str(actual), "--csv", "out.csv"]
        actual.unlink()
