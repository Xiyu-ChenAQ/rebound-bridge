"""Plot the saved fixed-timestep runtime-scaling measurements."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


PAPER_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path,
                        default=PAPER_ROOT / "results/scaling/complex_speed_ratio.csv")
    parser.add_argument("--output", type=Path,
                        default=PAPER_ROOT / "results/scaling/runtime_scaling_column")
    args = parser.parse_args()
    with args.csv.open(encoding="utf-8", newline="") as stream:
        rows = sorted(csv.DictReader(stream), key=lambda row: int(row["n_outer"]))
    if not rows:
        raise ValueError("No runtime measurements found")
    counts = [int(row["n_outer"]) for row in rows]
    ratios = []
    for row in rows:
        whfast = float(row["whfast_mean_seconds"])
        bridge = float(row["bridge_mean_seconds"])
        recorded_ratio = float(row["mean_speedup_whfast_over_bridge"])
        if not all(math.isfinite(value) and value > 0
                   for value in (whfast, bridge, recorded_ratio)):
            raise ValueError("Mean runtimes and recorded speedup must be finite and positive")
        ratio = whfast / bridge
        if not math.isclose(ratio, recorded_ratio, rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError("Recorded speedup disagrees with ratio of mean times")
        if (float(row["years"]) != 100 or float(row["dt_outer_days"]) != 4
                or row["bridge_method"] != "eos-lf4-2" or int(row["binary_ratio"]) != 64
                or float(row["whfast_dt_days"]) != 0.0625
                or int(row["repeats"]) != 3 or int(row["warmups"]) != 1):
            raise ValueError("Measurements do not share the stated configuration")
        ratios.append(ratio)

    plt.style.use(PAPER_ROOT / "mnras.mplstyle")
    fig, ax = plt.subplots(figsize=(3.35, 2.55), layout="constrained")
    colour = "#245c73"
    ax.plot(counts, ratios, color=colour, marker="o", label="Bridge LF(4,2)", zorder=3)
    ax.axhline(1, color="#777777", linestyle="--", linewidth=0.85, zorder=1)
    ax.set_xlim(-4, 136)
    ax.set_ylim(0, max(21.5, 1.15 * max(ratios)))
    ax.set_xticks([0, 32, 64, 96, 128])
    ax.set_yticks([0, 5, 10, 15, 20])
    ax.set_xlabel("Number of added bodies")
    ax.set_ylabel(r"Speedup $t_{\rm WHFast}/t_{\rm Bridge}$")
    ax.legend(loc="upper left")
    ax.annotate(rf"${ratios[0]:.2f}\times$", (counts[0], ratios[0]), xytext=(15, 5),
                textcoords="offset points", va="center", ha="left", color=colour,
                fontsize=8)
    ax.annotate(rf"${ratios[-1]:.1f}\times$", (counts[-1], ratios[-1]), xytext=(-4, 8),
                textcoords="offset points", va="bottom", ha="right", color=colour,
                fontsize=8.5)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output.with_suffix(".pdf"))
    fig.savefig(args.output.with_suffix(".png"))
    plt.close(fig)
    print(f"Saved {len(rows)} measured points; endpoint speedups {ratios[0]:.6f}, {ratios[-1]:.6f}")


if __name__ == "__main__":
    main()
