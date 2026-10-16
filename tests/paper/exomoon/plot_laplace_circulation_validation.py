"""Plot and summarize the validated loss of Laplace phase coherence."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "resonance_validation"
STYLE = ROOT / "mnras.mplstyle"

TIMESTEP_CASES = [
    ("Bridge 0.5 d", RESULTS / "a0p55_phase0_20yr_bridge_dt0p5.csv", "#b97820", "-"),
    ("Bridge 0.25 d", RESULTS / "a0p55_phase0_20yr_bridge_dt0p25.csv", "#36835b", "-"),
    ("Bridge 0.125 d", RESULTS / "a0p55_phase0_20yr_bridge_dt0p125.csv", "#245c73", "-"),
    ("IAS15", RESULTS / "a0p55_phase0_20yr_ias15.csv", "#333333", "--"),
]

PHASE_CASES = [
    (r"phase $0^\circ$", RESULTS / "a0p55_phase0_20yr_bridge_dt0p25.csv", "#245c73"),
    (r"phase $120^\circ$", RESULTS / "a0p55_phase120_20yr_bridge_dt0p25.csv", "#b97820"),
    (r"phase $240^\circ$", RESULTS / "a0p55_phase240_20yr_bridge_dt0p25.csv", "#777777"),
]


def read_laplace(path: Path) -> tuple[np.ndarray, np.ndarray, bool]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"empty CSV: {path}")
    time = np.asarray([float(row["time_yr"]) for row in rows])
    wrapped = np.asarray([float(row["laplace_deg"]) for row in rows])
    unwrapped = np.rad2deg(np.unwrap(np.deg2rad(wrapped)))
    escaped = any(float(row.get("ever_escaped", row.get("escaped", 0.0))) > 0.5 for row in rows)
    return time, unwrapped - unwrapped[0], escaped


def diagnostics(path: Path) -> dict[str, float]:
    time, angle, escaped = read_laplace(path)
    running_span = np.maximum.accumulate(angle) - np.minimum.accumulate(angle)
    crossings = np.flatnonzero(running_span >= 360.0)
    return {
        "final_time_yr": float(time[-1]),
        "laplace_span_deg": float(np.ptp(angle)),
        "laplace_net_turns": float(angle[-1] / 360.0),
        "max_sample_step_deg": float(np.max(np.abs(np.diff(angle)))),
        "first_full_circulation_yr": float(time[crossings[0]]) if crossings.size else float("nan"),
        "escaped": float(escaped),
    }


def write_summary(path: Path, data_dir: Path = RESULTS) -> None:
    rows = []
    for label, source, _color, _linestyle in TIMESTEP_CASES:
        rows.append({"suite": "integrator", "label": label, **diagnostics(data_dir / source.name)})
    for label, source, _color in PHASE_CASES[1:]:
        rows.append({"suite": "phase", "label": label, **diagnostics(data_dir / source.name)})
    fields = list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot(path: Path, data_dir: Path = RESULTS, *, column: bool = False) -> None:
    plt.style.use(STYLE)
    if column:
        fig, axes = plt.subplots(2, 1, figsize=(3.48, 5.4), constrained_layout=True)
    else:
        fig, axes = plt.subplots(1, 2, figsize=(7.03, 3.0), constrained_layout=True)

    for label, source, color, linestyle in TIMESTEP_CASES:
        time, angle, _escaped = read_laplace(data_dir / source.name)
        mask = time <= 5.0
        axes[0].plot(time[mask], angle[mask], color=color, linestyle=linestyle, label=label)
    ias15_crossing = diagnostics(data_dir / TIMESTEP_CASES[-1][1].name)["first_full_circulation_yr"]
    if np.isfinite(ias15_crossing):
        axes[0].axvline(ias15_crossing, color="0.55", linestyle=":", linewidth=0.8)
    axes[0].set_xlim(0.0, 5.0)
    axes[0].set_title("Integrator and timestep check")

    for label, source, color in PHASE_CASES:
        time, angle, _escaped = read_laplace(data_dir / source.name)
        axes[1].plot(time, angle, color=color, label=label)
    axes[1].set_xlim(0.0, 20.0)
    axes[1].set_title("Initial-phase dependence")

    for axis in axes:
        axis.axhline(-360.0, color="0.65", linestyle=":", linewidth=0.8)
        axis.set_xlabel("time [yr]")
        axis.set_ylabel(r"unwrapped $\phi_{\rm L}(t)-\phi_{\rm L}(0)$ [deg]")
        axis.legend(frameon=False, fontsize=7)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=RESULTS)
    parser.add_argument("--output-dir", type=Path, default=RESULTS)
    args = parser.parse_args()
    summary = args.output_dir / "laplace_circulation_validation_summary.csv"
    figure = args.output_dir / "laplace_circulation_validation.png"
    column_figure = figure.with_name(figure.stem + "_column.png")
    write_summary(summary, args.data_dir)
    plot(figure, args.data_dir)
    plot(column_figure, args.data_dir, column=True)
    print(f"wrote {summary}")
    print(f"wrote {figure}")
    print(f"wrote {figure.with_suffix('.pdf')}")
    print(f"wrote {column_figure}")
    print(f"wrote {column_figure.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
