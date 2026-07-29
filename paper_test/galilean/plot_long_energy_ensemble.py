"""Plot the four-realization cost-matched long-term ensemble."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "long_energy"
STYLE = ROOT / "mnras.mplstyle"
CASES = [
    {
        "label": "Bridge LF(4,2)",
        "path": BASE / "ensemble_lf42_bridge_10000yr_r4.csv",
        "color": "#245c73",
    },
    {
        "label": "corrected WHFast",
        "path": BASE / "ensemble_whfast_corr11c2_dt0p0085_10000yr_r4.csv",
        "color": "#8a4f18",
    },
]


def read_rows(path: Path) -> list[dict[str, float | str]]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"empty CSV: {path}")
    return [
        {key: value if key == "method" else float(value) for key, value in row.items()}
        for row in rows
    ]


def group_realizations(rows: list[dict[str, float | str]]) -> list[list[dict[str, float | str]]]:
    grouped: dict[int, list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        grouped[int(float(row["realization"]))].append(row)
    return [grouped[key] for key in sorted(grouped)]


def write_summary(path: Path) -> None:
    fields = [
        "label",
        "method",
        "realization",
        "final_time_yr",
        "final_wall_seconds",
        "max_energy_rel_abs",
        "final_energy_rel_abs",
        "max_angular_momentum_rel_abs",
        "max_linear_momentum_norm",
        "max_com_distance_au",
        "max_eccentricity",
        "escaped",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for case in CASES:
            for realization, rows in enumerate(group_realizations(read_rows(case["path"]))):
                writer.writerow(
                    {
                        "label": case["label"],
                        "method": rows[-1]["method"],
                        "realization": realization,
                        "final_time_yr": rows[-1]["time_yr"],
                        "final_wall_seconds": rows[-1]["wall_seconds"],
                        "max_energy_rel_abs": max(float(row["energy_rel_abs"]) for row in rows),
                        "final_energy_rel_abs": rows[-1]["energy_rel_abs"],
                        "max_angular_momentum_rel_abs": max(
                            float(row["angular_momentum_rel_abs"]) for row in rows
                        ),
                        "max_linear_momentum_norm": max(float(row["linear_momentum_norm"]) for row in rows),
                        "max_com_distance_au": max(float(row["com_distance_au"]) for row in rows),
                        "max_eccentricity": max(float(row["max_eccentricity"]) for row in rows),
                        "escaped": int(any(float(row["escaped"]) > 0.5 for row in rows)),
                    }
                )


def plot(path: Path) -> None:
    plt.style.use(STYLE)
    fig, axes = plt.subplots(1, 2, figsize=(7.03, 3.0), constrained_layout=True)
    for case in CASES:
        realizations = group_realizations(read_rows(case["path"]))
        for rows in realizations:
            axes[0].plot(
                [float(row["time_yr"]) for row in rows],
                [float(row["energy_rel_abs"]) for row in rows],
                color=case["color"],
                linewidth=0.55,
                alpha=0.32,
            )
            axes[1].plot(
                [float(row["time_yr"]) for row in rows],
                [float(row["angular_momentum_rel_abs"]) for row in rows],
                color=case["color"],
                linewidth=0.55,
                alpha=0.32,
            )

        times = np.asarray([float(row["time_yr"]) for row in realizations[0]])
        energy = np.asarray(
            [[float(row["energy_rel_abs"]) for row in rows] for rows in realizations]
        )
        angular = np.asarray(
            [[float(row["angular_momentum_rel_abs"]) for row in rows] for rows in realizations]
        )
        axes[0].plot(times, np.mean(energy, axis=0), color=case["color"], linewidth=1.8, label=case["label"])
        axes[1].plot(times, np.mean(angular, axis=0), color=case["color"], linewidth=1.8, label=case["label"])

    for ax in axes:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("time [yr]")
    axes[0].set_ylabel("|relative energy error|")
    axes[1].set_ylabel("|relative angular-momentum error|")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    summary = BASE / "galilean_long_ensemble_10000yr_r4_summary.csv"
    figure = BASE / "galilean_long_ensemble_10000yr_r4.png"
    write_summary(summary)
    plot(figure)
    print(f"wrote {summary}")
    print(f"wrote {figure}")
    print(f"wrote {figure.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
