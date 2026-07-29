"""Plot and summarize the LF(4,2) generalized-order crossover."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "generalized_order"
DEFAULT_INPUT = BASE / "generalized_order_crossover.csv"
STYLE = ROOT / "mnras.mplstyle"


def resolved_mask(errors: list[float], reference_floor: float) -> list[bool]:
    return [error > 10.0 * reference_floor for error in errors]


def local_orders(steps: list[float], errors: list[float]) -> tuple[list[float], list[float]]:
    midpoints: list[float] = []
    orders: list[float] = []
    for h0, h1, e0, e1 in zip(steps, steps[1:], errors, errors[1:]):
        midpoints.append(math.sqrt(h0 * h1))
        orders.append(math.log(e1 / e0) / math.log(h1 / h0))
    return midpoints, orders


def read_rows(path: Path) -> list[dict[str, float | str]]:
    with path.open(newline="") as stream:
        raw = list(csv.DictReader(stream))
    if not raw:
        raise RuntimeError(f"empty crossover table: {path}")
    rows: list[dict[str, float | str]] = []
    for row in raw:
        rows.append(
            {
                key: value if key == "method" else float(value)
                for key, value in row.items()
            }
        )
    return rows


def grouped_rows(
    rows: list[dict[str, float | str]],
) -> dict[tuple[str, float], list[dict[str, float | str]]]:
    groups: dict[tuple[str, float], list[dict[str, float | str]]] = defaultdict(list)
    for row in rows:
        groups[(str(row["method"]), float(row["coupling_scale"]))].append(row)
    for group in groups.values():
        group.sort(key=lambda row: float(row["dt_days"]))
    return groups


def fit_order(steps: list[float], errors: list[float]) -> float:
    return float(np.polyfit(np.log(steps), np.log(errors), 1)[0])


def write_summary(rows: list[dict[str, float | str]], path: Path) -> None:
    fields = [
        "coupling_scale",
        "effective_epsilon",
        "method",
        "resolved_points",
        "fit_order",
        "fine_local_order",
        "coarse_local_order",
        "reference_floor_au",
    ]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for (method, coupling_scale), group in sorted(grouped_rows(rows).items(), key=lambda item: (item[0][1], item[0][0])):
            steps = [float(row["dt_days"]) for row in group]
            errors = [float(row["max_position_error_au"]) for row in group]
            floor = float(group[0]["reference_floor_au"])
            mask = resolved_mask(errors, floor)
            resolved_steps = [step for step, keep in zip(steps, mask) if keep]
            resolved_errors = [error for error, keep in zip(errors, mask) if keep]
            local_midpoints, orders = local_orders(steps, errors)
            local_resolved = [
                order
                for index, order in enumerate(orders)
                if mask[index] and mask[index + 1]
            ]
            writer.writerow(
                {
                    "coupling_scale": coupling_scale,
                    "effective_epsilon": group[0]["effective_epsilon"],
                    "method": method,
                    "resolved_points": len(resolved_steps),
                    "fit_order": fit_order(resolved_steps, resolved_errors),
                    "fine_local_order": local_resolved[0] if local_resolved else "",
                    "coarse_local_order": local_resolved[-1] if local_resolved else "",
                    "reference_floor_au": floor,
                }
            )


def plot(rows: list[dict[str, float | str]], png_path: Path, pdf_path: Path) -> None:
    groups = grouped_rows(rows)
    scales = [0.01, 0.1, 1.0, 5.0]
    local_order_scales = {0.01, 5.0}
    colours = ["#173f73", "#5f6b7a", "#8c8061", "#c46b28"]

    plt.style.use(STYLE)
    figure, axes = plt.subplots(1, 2, figsize=(7.03, 2.85), constrained_layout=True)

    for coupling_scale, colour in zip(scales, colours):
        group = groups[("lf42", coupling_scale)]
        steps = [float(row["dt_days"]) for row in group]
        errors = [float(row["max_position_error_au"]) for row in group]
        epsilon = float(group[0]["effective_epsilon"])
        floor = float(group[0]["reference_floor_au"])
        mask = resolved_mask(errors, floor)
        exponent = math.floor(math.log10(epsilon))
        mantissa = epsilon / 10.0**exponent
        label = rf"$\epsilon={mantissa:.1f}\times10^{{{exponent}}}$"
        axes[0].plot(steps, errors, marker="o", color=colour, label=label)

        if coupling_scale in local_order_scales:
            midpoints, orders = local_orders(steps, errors)
            pair_mask = [mask[index] and mask[index + 1] for index in range(len(mask) - 1)]
            axes[1].plot(
                [midpoint for midpoint, keep in zip(midpoints, pair_mask) if keep],
                [order for order, keep in zip(orders, pair_mask) if keep],
                marker="o",
                color=colour,
                label=label,
            )

    axes[0].set_xscale("log", base=2)
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"outer step $h$ [d]")
    axes[0].set_ylabel("maximum position error [AU]")
    axes[0].set_title("Error")

    axes[1].set_xscale("log", base=2)
    axes[1].axhline(2.0, color="0.45", linestyle="--", linewidth=1.0, label="order 2")
    axes[1].axhline(4.0, color="0.45", linestyle=":", linewidth=1.0, label="order 4")
    axes[1].set_xlabel(r"geometric-mean step $h$ [d]")
    axes[1].set_ylabel("local convergence order")
    axes[1].set_ylim(1.0, 5.0)
    axes[1].set_title("Local order")

    axes[0].legend(ncol=2)
    axes[1].legend(ncol=2)
    figure.savefig(png_path)
    figure.savefig(pdf_path)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()

    rows = read_rows(args.input)
    summary = BASE / "generalized_order_crossover_summary.csv"
    png = BASE / "generalized_order_crossover.png"
    pdf = BASE / "generalized_order_crossover.pdf"
    write_summary(rows, summary)
    plot(rows, png, pdf)
    print(f"wrote {summary}")
    print(f"wrote {png}")
    print(f"wrote {pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
