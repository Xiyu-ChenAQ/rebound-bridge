"""Run and plot Bridge outer/inner timestep convergence separately."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


AU_KM = 149597870.7
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "two_timestep_convergence"
DEFAULT_EXE = ROOT / "build" / "Release" / "galilean_horizons_compare.exe"
METHODS = [("order2", "KDK / order2"), ("eos_lf4_2", "EOS LF4_2")]
STYLE = ROOT / "mnras.mplstyle"
OUTER_DTS = [0.125, 0.25, 0.5, 1.0]
INNER_RATIOS = [5, 10, 20, 40, 80, 160, 320]
YEARS = 4.0


def token(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def case_path(sweep: str, method: str, dt_outer_days: float, inner_ratio: int) -> Path:
    return RESULTS / (
        f"galilean_{sweep}_{method}_{token(YEARS)}yr_"
        f"dt{token(dt_outer_days)}_ir{inner_ratio}.csv"
    )


def run_case(exe: Path, sweep: str, method: str, dt_outer_days: float, inner_ratio: int) -> Path:
    path = case_path(sweep, method, dt_outer_days, inner_ratio)
    path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(exe),
        "--method",
        method,
        "--years",
        f"{YEARS:.17g}",
        "--samples",
        "1",
        "--dt-outer-days",
        f"{dt_outer_days:.17g}",
        "--inner-ratio",
        str(inner_ratio),
        "--csv",
        str(path),
    ]
    print("running", " ".join(command))
    subprocess.run(command, check=True)
    return path


def read_last(path: Path) -> dict[str, float]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"empty CSV: {path}")
    return {key: float(value) for key, value in rows[-1].items()}


def collect_rows(exe: Path, skip_run: bool) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    cases: list[tuple[str, str, str, float, int]] = []
    for method, label in METHODS:
        for dt_outer in OUTER_DTS:
            cases.append(("outer", method, label, dt_outer, 200))
        for inner_ratio in INNER_RATIOS:
            cases.append(("inner", method, label, 0.5, inner_ratio))

    for sweep, method, label, dt_outer, inner_ratio in cases:
        path = case_path(sweep, method, dt_outer, inner_ratio)
        if not skip_run:
            path = run_case(exe, sweep, method, dt_outer, inner_ratio)
        last = read_last(path)
        rows.append(
            {
                "sweep": sweep,
                "method": label,
                "dt_outer_days": dt_outer,
                "inner_ratio": float(inner_ratio),
                "dt_inner_days": dt_outer / inner_ratio,
                "moon_error_km": last["max_moon_error_au"] * AU_KM,
                "laplace_error_deg": last["laplace_diff_deg"],
                "energy_rel_abs": abs(last["bridge_energy_rel"]),
                "speedup": last["cumulative_speedup"],
            }
        )
    return rows


def fit_outer_slope(rows: list[dict[str, float | str]], method: str, key: str) -> float:
    selected = [
        row
        for row in rows
        if row["sweep"] == "outer"
        and row["method"] == method
        and float(row["dt_outer_days"]) <= 0.5
    ]
    x = np.log10([float(row["dt_outer_days"]) for row in selected])
    y = np.log10([float(row[key]) for row in selected])
    return float(np.polyfit(x, y, 1)[0])


def write_summary(rows: list[dict[str, float | str]], path: Path) -> None:
    fields = [
        "sweep",
        "method",
        "dt_outer_days",
        "inner_ratio",
        "dt_inner_days",
        "moon_error_km",
        "laplace_error_deg",
        "energy_rel_abs",
        "speedup",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot(rows: list[dict[str, float | str]], path: Path) -> None:
    styles = {
        "KDK / order2": {"marker": "s", "color": "#a03c30", "label": "KDK"},
        "EOS LF4_2": {"marker": "o", "color": "#245c73", "label": "LF(4,2)"},
    }
    plt.style.use(STYLE)
    fig, axes = plt.subplots(2, 2, figsize=(7.03, 5.05), constrained_layout=True)

    for method, style in styles.items():
        outer = sorted(
            [row for row in rows if row["sweep"] == "outer" and row["method"] == method],
            key=lambda row: float(row["dt_outer_days"]),
        )
        inner = sorted(
            [row for row in rows if row["sweep"] == "inner" and row["method"] == method],
            key=lambda row: float(row["dt_inner_days"]),
        )
        outer_x = [float(row["dt_outer_days"]) for row in outer]
        inner_x = [float(row["dt_inner_days"]) for row in inner]
        moon_slope = fit_outer_slope(rows, method, "moon_error_km")
        laplace_slope = fit_outer_slope(rows, method, "laplace_error_deg")

        axes[0, 0].plot(
            outer_x,
            [float(row["moon_error_km"]) for row in outer],
            label=f"{style['label']} ($p={moon_slope:.2f}$)",
            marker=style["marker"],
            color=style["color"],
        )
        axes[0, 1].plot(
            outer_x,
            [float(row["laplace_error_deg"]) for row in outer],
            label=f"{style['label']} ($p={laplace_slope:.2f}$)",
            marker=style["marker"],
            color=style["color"],
        )
        axes[1, 0].plot(
            inner_x,
            [float(row["moon_error_km"]) for row in inner],
            label=style["label"],
            marker=style["marker"],
            color=style["color"],
        )
        axes[1, 1].plot(
            inner_x,
            [float(row["laplace_error_deg"]) for row in inner],
            label=style["label"],
            marker=style["marker"],
            color=style["color"],
        )

    for ax in axes.ravel():
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.legend()

    axes[0, 0].set_title("Outer step: position")
    axes[0, 1].set_title("Outer step: phase")
    axes[0, 0].set_ylabel("position error [km]")
    axes[0, 1].set_ylabel("Laplace-angle error [deg]")
    axes[0, 0].set_xlabel(r"outer step $h$ [d]")
    axes[0, 1].set_xlabel(r"outer step $h$ [d]")

    axes[1, 0].set_title("Inner step: position")
    axes[1, 1].set_title("Inner step: phase")
    axes[1, 0].set_ylabel("position error [km]")
    axes[1, 1].set_ylabel("Laplace-angle error [deg]")
    axes[1, 0].set_xlabel(r"inner step $\delta$ [d]")
    axes[1, 1].set_xlabel(r"inner step $\delta$ [d]")

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--skip-run", action="store_true")
    args = parser.parse_args()

    rows = collect_rows(args.exe, args.skip_run)
    summary = RESULTS / "galilean_two_timestep_convergence_summary.csv"
    figure = RESULTS / "galilean_two_timestep_convergence.png"
    write_summary(rows, summary)
    plot(rows, figure)
    print(f"wrote {summary}")
    print(f"wrote {figure}")
    print(f"wrote {figure.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
