from __future__ import annotations

from pathlib import Path
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
