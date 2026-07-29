"""Run and aggregate repeated KDK/LF(4,2) Galilean benchmarks."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import math
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import time


LF42_A = 0.2113248654051871
INNER_DT_DAYS = 0.0125
OUTER_DT_DAYS = (0.125, 0.25, 0.5, 1.0)


@dataclass(frozen=True)
class Case:
    method: str
    dt_outer_days: float
    inner_ratio: int
    dt_inner_days: float


@dataclass(frozen=True)
class WorkCounts:
    outer_steps: int
    flow_segments: int
    cross_kicks: int
    internal_steps: int


@dataclass(frozen=True)
class ScheduledCase:
    repeat: int
    order: int
    case: Case


@dataclass(frozen=True)
class RunMetrics:
    method_seconds: float
    ias15_seconds: float
    speedup: float
    max_moon_error_au: float
    max_laplace_diff_deg: float


@dataclass(frozen=True)
class AggregateSummary:
    case: Case
    repeats: int
    method_seconds_median: float
    method_seconds_q25: float
    method_seconds_q75: float
    ias15_seconds_median: float
    ias15_seconds_q25: float
    ias15_seconds_q75: float
    speedup_median: float
    speedup_q25: float
    speedup_q75: float
    max_moon_error_au: float
    max_laplace_diff_deg: float


def _exact_integer(value: float, name: str) -> int:
    rounded = round(value)
    if not math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-10):
        raise ValueError(f"{name} must be an integer, got {value:.17g}")
    return int(rounded)


def build_cases() -> list[Case]:
    cases: list[Case] = []
    for method in ("order2", "eos_lf4_2"):
        for dt_outer_days in OUTER_DT_DAYS:
            inner_ratio = _exact_integer(
                dt_outer_days / INNER_DT_DAYS,
                "outer-to-inner ratio",
            )
            cases.append(Case(method, dt_outer_days, inner_ratio, INNER_DT_DAYS))
    return cases


def build_run_schedule(repeats: int, seed: int) -> list[ScheduledCase]:
    if repeats < 1:
        raise ValueError("repeats must be positive")
    schedule: list[ScheduledCase] = []
    base_cases = build_cases()
    for repeat in range(1, repeats + 1):
        round_cases = list(base_cases)
        random.Random(seed + repeat).shuffle(round_cases)
        for order, case in enumerate(round_cases, start=1):
            schedule.append(ScheduledCase(repeat=repeat, order=order, case=case))
    return schedule


def validate_configuration(years: float, output_interval_days: float, cases: list[Case]) -> int:
    if years <= 0.0:
        raise ValueError("years must be positive")
    if output_interval_days <= 0.0:
        raise ValueError("output interval must be positive")
    duration_days = years * 365.25
    output_intervals = _exact_integer(
        duration_days / output_interval_days,
        "output intervals",
    )
    for case in cases:
        _exact_integer(
            output_interval_days / case.dt_outer_days,
            "output-to-outer ratio",
        )
        _exact_integer(
            case.dt_outer_days / case.dt_inner_days,
            "outer-to-inner ratio",
        )
    return output_intervals


def compute_work_counts(
    method: str,
    duration_days: float,
    output_intervals: int,
    dt_outer_days: float,
    dt_inner_days: float,
) -> WorkCounts:
    outer_steps = _exact_integer(duration_days / dt_outer_days, "outer steps")
    inner_ratio = _exact_integer(dt_outer_days / dt_inner_days, "inner ratio")

    if method == "order2":
        return WorkCounts(
            outer_steps=outer_steps,
            flow_segments=outer_steps,
            cross_kicks=outer_steps + output_intervals,
            internal_steps=outer_steps * inner_ratio,
        )
    if method == "eos_lf4_2":
        segment_steps = (
            2 * math.ceil(LF42_A * inner_ratio)
            + math.ceil((1.0 - 2.0 * LF42_A) * inner_ratio)
        )
        return WorkCounts(
            outer_steps=outer_steps,
            flow_segments=3 * outer_steps,
            cross_kicks=2 * outer_steps,
            internal_steps=outer_steps * segment_steps,
        )
    raise ValueError(f"unsupported Bridge method: {method}")


def inclusive_quantile(values: list[float], quantile: float) -> float:
    if not values:
        raise ValueError("cannot calculate a quantile of an empty sequence")
    if not 0.0 <= quantile <= 1.0:
        raise ValueError("quantile must be between zero and one")

    ordered = sorted(values)
    index = (len(ordered) - 1) * quantile
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    fraction = index - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def _read_run_metrics(path: Path, expected_years: float, expected_samples: int) -> RunMetrics:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != expected_samples:
        raise ValueError(f"wrong sample count in {path}: {len(rows)}")

    final = rows[-1]
    final_time = float(final["time_yr"])
    if not math.isclose(final_time, expected_years, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError(f"wrong final time in {path}: {final_time:.17g}")
    if any(row["bridge_escaped"] == "1" or row["ias15_escaped"] == "1" for row in rows):
        raise ValueError(f"escaped run in {path}")

    method_seconds = float(final["bridge_cumulative_seconds"])
    ias15_seconds = float(final["ias15_cumulative_seconds"])
    speedup = ias15_seconds / method_seconds
    values = (method_seconds, ias15_seconds, speedup)
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError(f"invalid timing in {path}")

    return RunMetrics(
        method_seconds=method_seconds,
        ias15_seconds=ias15_seconds,
        speedup=speedup,
        max_moon_error_au=max(float(row["max_moon_error_au"]) for row in rows),
        max_laplace_diff_deg=max(float(row["laplace_diff_deg"]) for row in rows),
    )


def _require_repeat_stability(values: list[float], label: str) -> float:
    reference = values[0]
    for value in values[1:]:
        if not math.isclose(value, reference, rel_tol=1.0e-13, abs_tol=1.0e-15):
            raise ValueError(f"repeat-dependent {label}: {values}")
    return reference


def aggregate_case(
    case: Case,
    paths: list[Path],
    expected_years: float,
    expected_samples: int,
) -> AggregateSummary:
    if not paths:
        raise ValueError("cannot aggregate an empty case")
    runs = [_read_run_metrics(path, expected_years, expected_samples) for path in paths]
    method_times = [run.method_seconds for run in runs]
    ias15_times = [run.ias15_seconds for run in runs]
    speedups = [run.speedup for run in runs]
    moon_errors = [run.max_moon_error_au for run in runs]
    laplace_errors = [run.max_laplace_diff_deg for run in runs]

    return AggregateSummary(
        case=case,
        repeats=len(runs),
        method_seconds_median=statistics.median(method_times),
        method_seconds_q25=inclusive_quantile(method_times, 0.25),
        method_seconds_q75=inclusive_quantile(method_times, 0.75),
        ias15_seconds_median=statistics.median(ias15_times),
        ias15_seconds_q25=inclusive_quantile(ias15_times, 0.25),
        ias15_seconds_q75=inclusive_quantile(ias15_times, 0.75),
        speedup_median=statistics.median(speedups),
        speedup_q25=inclusive_quantile(speedups, 0.25),
        speedup_q75=inclusive_quantile(speedups, 0.75),
        max_moon_error_au=_require_repeat_stability(moon_errors, "moon error"),
        max_laplace_diff_deg=_require_repeat_stability(laplace_errors, "Laplace error"),
    )


def _method_label(method: str) -> str:
    if method == "order2":
        return "KDK / order2"
    if method == "eos_lf4_2":
        return "EOS LF4_2"
    raise ValueError(f"unsupported Bridge method: {method}")


def _number_slug(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def _raw_paths(output_dir: Path, item: ScheduledCase) -> tuple[Path, Path, Path]:
    stem = (
        f"repeat{item.repeat:02d}_order{item.order:02d}_{item.case.method}_"
        f"dt{_number_slug(item.case.dt_outer_days)}_ir{item.case.inner_ratio}"
    )
    raw_dir = output_dir / "raw"
    return (
        raw_dir / f"{stem}.csv",
        raw_dir / f"{stem}.stdout.txt",
        raw_dir / f"{stem}.stderr.txt",
    )


def _runner_command(
    executable: Path,
    case: Case,
    years: float,
    samples: int,
    csv_path: Path,
) -> list[str]:
    return [
        str(executable),
        "--method",
        case.method,
        "--years",
        f"{years:.17g}",
        "--samples",
        str(samples),
        "--dt-outer-days",
        f"{case.dt_outer_days:.17g}",
        "--inner-ratio",
        str(case.inner_ratio),
        "--csv",
        str(csv_path),
    ]


def _run_command(command: list[str], stdout_path: Path, stderr_path: Path) -> float:
    started = time.perf_counter()
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    elapsed = time.perf_counter() - started
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"benchmark command failed with exit code {completed.returncode}: {' '.join(command)}"
        )
    return elapsed


def _run_warmups(executable: Path, output_dir: Path) -> None:
    warmup_dir = output_dir / "warmup"
    warmup_dir.mkdir(parents=True, exist_ok=True)
    for method in ("order2", "eos_lf4_2"):
        case = next(
            case
            for case in build_cases()
            if case.method == method and case.dt_outer_days == 0.5
        )
        csv_path = warmup_dir / f"{method}.csv"
        stdout_path = warmup_dir / f"{method}.stdout.txt"
        stderr_path = warmup_dir / f"{method}.stderr.txt"
        command = _runner_command(
            executable,
            case,
            years=5.0 / 365.25,
            samples=1,
            csv_path=csv_path,
        )
        _run_command(command, stdout_path, stderr_path)


def _summary_fields() -> list[str]:
    return [
        "method",
        "method_label",
        "dt_outer_days",
        "inner_ratio",
        "dt_inner_days",
        "repeats",
        "method_seconds_median",
        "method_seconds_q25",
        "method_seconds_q75",
        "ias15_seconds_median",
        "ias15_seconds_q25",
        "ias15_seconds_q75",
        "speedup_median",
        "speedup_q25",
        "speedup_q75",
        "max_moon_error_au",
        "max_laplace_diff_deg",
    ]


def _write_summary(path: Path, summaries: list[AggregateSummary]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=_summary_fields())
        writer.writeheader()
        for summary in summaries:
            case = summary.case
            writer.writerow(
                {
                    "method": case.method,
                    "method_label": _method_label(case.method),
                    "dt_outer_days": case.dt_outer_days,
                    "inner_ratio": case.inner_ratio,
                    "dt_inner_days": case.dt_inner_days,
                    "repeats": summary.repeats,
                    "method_seconds_median": summary.method_seconds_median,
                    "method_seconds_q25": summary.method_seconds_q25,
                    "method_seconds_q75": summary.method_seconds_q75,
                    "ias15_seconds_median": summary.ias15_seconds_median,
                    "ias15_seconds_q25": summary.ias15_seconds_q25,
                    "ias15_seconds_q75": summary.ias15_seconds_q75,
                    "speedup_median": summary.speedup_median,
                    "speedup_q25": summary.speedup_q25,
                    "speedup_q75": summary.speedup_q75,
                    "max_moon_error_au": summary.max_moon_error_au,
                    "max_laplace_diff_deg": summary.max_laplace_diff_deg,
                }
            )


def _write_work_counts(
    path: Path,
    cases: list[Case],
    duration_days: float,
    output_intervals: int,
) -> None:
    fields = [
        "method",
        "method_label",
        "dt_outer_days",
        "inner_ratio",
        "dt_inner_days",
        "duration_days",
        "output_intervals",
        "outer_steps",
        "flow_segments",
        "cross_kicks",
        "internal_steps",
    ]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for case in cases:
            counts = compute_work_counts(
                case.method,
                duration_days,
                output_intervals,
                case.dt_outer_days,
                case.dt_inner_days,
            )
            writer.writerow(
                {
                    "method": case.method,
                    "method_label": _method_label(case.method),
                    "dt_outer_days": case.dt_outer_days,
                    "inner_ratio": case.inner_ratio,
                    "dt_inner_days": case.dt_inner_days,
                    "duration_days": duration_days,
                    "output_intervals": output_intervals,
                    "outer_steps": counts.outer_steps,
                    "flow_segments": counts.flow_segments,
                    "cross_kicks": counts.cross_kicks,
                    "internal_steps": counts.internal_steps,
                }
            )


def _read_manifest_metrics(
    item: ScheduledCase,
    csv_path: Path,
    years: float,
    samples: int,
) -> dict[str, object]:
    metrics = _read_run_metrics(csv_path, years, samples)
    return {
        "repeat": item.repeat,
        "order": item.order,
        "method": item.case.method,
        "dt_outer_days": item.case.dt_outer_days,
        "inner_ratio": item.case.inner_ratio,
        "dt_inner_days": item.case.dt_inner_days,
        "raw_csv": str(csv_path.resolve()),
        "method_seconds": metrics.method_seconds,
        "ias15_seconds": metrics.ias15_seconds,
        "speedup": metrics.speedup,
        "max_moon_error_au": metrics.max_moon_error_au,
        "max_laplace_diff_deg": metrics.max_laplace_diff_deg,
    }


def _write_manifest(path: Path, records: list[dict[str, object]]) -> None:
    fields = [
        "repeat",
        "order",
        "started_utc",
        "method",
        "dt_outer_days",
        "inner_ratio",
        "dt_inner_days",
        "raw_csv",
        "stdout_log",
        "stderr_log",
        "process_elapsed_seconds",
        "method_seconds",
        "ias15_seconds",
        "speedup",
        "max_moon_error_au",
        "max_laplace_diff_deg",
        "executable",
        "platform",
        "python_version",
        "build_config",
    ]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def _aggregate_all(
    output_dir: Path,
    schedule: list[ScheduledCase],
    years: float,
    samples: int,
) -> list[AggregateSummary]:
    paths_by_case: dict[Case, list[Path]] = {case: [] for case in build_cases()}
    for item in schedule:
        csv_path, _stdout_path, _stderr_path = _raw_paths(output_dir, item)
        paths_by_case[item.case].append(csv_path)
    return [
        aggregate_case(case, paths_by_case[case], years, samples)
        for case in build_cases()
    ]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--years", type=float, default=20.0)
    parser.add_argument("--output-interval-days", type=float, default=5.0)
    parser.add_argument("--seed", type=int, default=20260713)
    parser.add_argument("--build-config", default="Release")
    parser.add_argument("--skip-run", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    cases = build_cases()
    samples = validate_configuration(args.years, args.output_interval_days, cases)
    duration_days = args.years * 365.25
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "raw").mkdir(parents=True, exist_ok=True)
    schedule = build_run_schedule(args.repeats, args.seed)

    executable = args.exe.resolve()
    if not args.skip_run:
        if not executable.is_file():
            raise FileNotFoundError(f"benchmark executable not found: {executable}")
        _run_warmups(executable, output_dir)
        manifest_records: list[dict[str, object]] = []
        for index, item in enumerate(schedule, start=1):
            csv_path, stdout_path, stderr_path = _raw_paths(output_dir, item)
            started_utc = datetime.now(timezone.utc).isoformat()
            command = _runner_command(executable, item.case, args.years, samples, csv_path)
            print(
                f"[{index:02d}/{len(schedule)}] repeat={item.repeat} order={item.order} "
                f"method={item.case.method} h={item.case.dt_outer_days:g} d",
                flush=True,
            )
            process_elapsed = _run_command(command, stdout_path, stderr_path)
            record = _read_manifest_metrics(item, csv_path, args.years, samples)
            record.update(
                {
                    "started_utc": started_utc,
                    "stdout_log": str(stdout_path.resolve()),
                    "stderr_log": str(stderr_path.resolve()),
                    "process_elapsed_seconds": process_elapsed,
                    "executable": str(executable),
                    "platform": platform.platform(),
                    "python_version": sys.version.split()[0],
                    "build_config": args.build_config,
                }
            )
            manifest_records.append(record)
        _write_manifest(output_dir / "run_manifest.csv", manifest_records)

    summaries = _aggregate_all(output_dir, schedule, args.years, samples)
    _write_summary(output_dir / "kdk_lf42_cost_accuracy_summary.csv", summaries)
    _write_work_counts(
        output_dir / "kdk_lf42_work_counts.csv",
        cases,
        duration_days,
        samples,
    )
    print(f"wrote {output_dir / 'kdk_lf42_cost_accuracy_summary.csv'}")
    print(f"wrote {output_dir / 'kdk_lf42_work_counts.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
