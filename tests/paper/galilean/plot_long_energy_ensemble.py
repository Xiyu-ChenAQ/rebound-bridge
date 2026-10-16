"""Plot the four-realization long-term conservation ensemble."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "long_energy"
STYLE = ROOT / "mnras.mplstyle"
REALIZATIONS = 4
FINAL_TIME_YR = 10000.0
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
        realization = float(row["realization"])
        if not np.isfinite(realization) or not realization.is_integer():
            raise ValueError(f"invalid realization ID: {realization}")
        grouped[int(realization)].append(row)
    if set(grouped) != set(range(REALIZATIONS)):
        raise ValueError(f"expected realization IDs 0 through {REALIZATIONS - 1}, got {sorted(grouped)}")

    reference_times = None
    for realization, samples in sorted(grouped.items()):
        times = np.asarray([float(row["time_yr"]) for row in samples])
        if len(times) < 2 or not np.all(np.isfinite(times)) or np.any(times <= 0) or np.any(np.diff(times) <= 0):
            raise ValueError(f"realization {realization}: times must be finite, positive and strictly increasing")
        if reference_times is not None and not np.array_equal(times, reference_times):
            raise ValueError(f"realization {realization}: output times do not match realization 0")
        reference_times = times
        if float(samples[0]["step"]) != 1:
            raise ValueError(f"realization {realization}: missing first integration step")
        # The benchmark rounds the requested duration to the nearest full step.
        dt_key = "dt_outer_days" if samples[0]["method"] == "bridge" else "dt_days"
        dt_years = float(samples[0][dt_key]) / 365.25
        if not np.isfinite(dt_years) or dt_years <= 0 or abs(times[-1] - FINAL_TIME_YR) > dt_years / 2 + 1e-10:
            raise ValueError(f"realization {realization}: incomplete {FINAL_TIME_YR:g}-year integration")
        for field in ("energy_rel_abs", "angular_momentum_rel_abs", "wall_seconds", "linear_momentum_norm", "com_distance_au", "max_eccentricity"):
            values = np.asarray([float(row[field]) for row in samples])
            if not np.all(np.isfinite(values)) or np.any(values < 0):
                raise ValueError(f"realization {realization}: {field} must be finite and nonnegative")
    return [grouped[key] for key in sorted(grouped)]


def load_cases(data_dir: Path) -> list[dict]:
    cases = []
    for case in CASES:
        path = data_dir / case["path"].name
        try:
            realizations = group_realizations(read_rows(path))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid ensemble CSV {path}: {exc}") from exc
        cases.append({**case, "realizations": realizations})
    return cases


def write_summary(path: Path, cases: list[dict]) -> None:
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
        for case in cases:
            for rows in case["realizations"]:
                writer.writerow(
                    {
                        "label": case["label"],
                        "method": rows[-1]["method"],
                        "realization": int(rows[-1]["realization"]),
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


def plot(path: Path, cases: list[dict], *, column: bool = False) -> None:
    plt.style.use(STYLE)
    fig, axes = plt.subplots(
        2 if column else 1,
        1 if column else 2,
        figsize=(3.35, 4.2) if column else (7.03, 3.0),
        sharex=column,
        constrained_layout=True,
    )
    for case in cases:
        realizations = case["realizations"]
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
        if not column or ax is axes[-1]:
            ax.set_xlabel("time [yr]")
    axes[0].set_ylabel(r"$|\Delta E|/|E_0|$")
    axes[1].set_ylabel(r"$\left||\mathbf{L}|-|\mathbf{L}_0|\right|/|\mathbf{L}_0|$")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=BASE, help="directory containing the two ensemble CSVs")
    parser.add_argument("--output-dir", type=Path, default=BASE, help="directory for the summary and figures")
    args = parser.parse_args()
    cases = load_cases(args.data_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = args.output_dir / "galilean_long_ensemble_10000yr_r4_summary.csv"
    write_summary(summary, cases)
    print(f"wrote {summary}")
    for column in (False, True):
        suffix = "_column" if column else ""
        figure = args.output_dir / f"galilean_long_ensemble_10000yr_r4{suffix}.png"
        plot(figure, cases, column=column)
        print(f"wrote {figure}")
        print(f"wrote {figure.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
