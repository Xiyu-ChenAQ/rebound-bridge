#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_rows(path):
    with path.open(newline="", encoding="utf-8") as fp:
        reader = csv.DictReader(fp)
        rows = []
        for row in reader:
            rows.append(
                {
                    "years": float(row["years"]),
                    "samples": int(row["samples"]),
                    "bridge_total_seconds": float(row["bridge_total_seconds"]),
                    "ias15_total_seconds": float(row["ias15_total_seconds"]),
                    "overall_speedup": float(row["overall_speedup"]),
                }
            )
    return rows


def save_fig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_sweep(rows, outdir):
    years = [row["years"] for row in rows]
    bridge_total = [row["bridge_total_seconds"] for row in rows]
    ias15_total = [row["ias15_total_seconds"] for row in rows]
    bridge_ms = [1e3 * row["bridge_total_seconds"] / row["samples"] for row in rows]
    ias15_ms = [1e3 * row["ias15_total_seconds"] / row["samples"] for row in rows]
    speedup = [row["overall_speedup"] for row in rows]

    fig, axes = plt.subplots(3, 1, figsize=(12, 13), sharex=True)
    ax = axes[0]
    ax.loglog(years, bridge_total, marker="o", label="Bridge")
    ax.loglog(years, ias15_total, marker="o", label="IAS15")
    ax.set_ylabel("Total runtime [s]")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    ax.semilogx(years, bridge_ms, marker="o", label="Bridge")
    ax.semilogx(years, ias15_ms, marker="o", label="IAS15")
    ax.set_ylabel("Mean time / sample [ms]")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[2]
    ax.semilogx(years, speedup, marker="o", color="tab:red")
    ax.axhline(1.0, color="0.3", lw=0.8)
    ax.set_ylabel("IAS15 / Bridge")
    ax.set_xlabel("Years")
    ax.grid(True, which="both", alpha=0.25)

    save_fig(fig, outdir / "efficiency_sweep.png")


def write_summary(rows, outdir):
    best = max(rows, key=lambda row: row["overall_speedup"])
    worst = min(rows, key=lambda row: row["overall_speedup"])
    longest = max(rows, key=lambda row: row["years"])

    lines = [
        "Bridge vs IAS15 years sweep summary",
        "speedup definition: IAS15 total runtime / Bridge total runtime; values > 1 mean Bridge is faster",
        f"cases: {len(rows)}",
        f"years_min: {rows[0]['years']:.3f}",
        f"years_max: {rows[-1]['years']:.3f}",
        f"best_speedup_years: {best['years']:.3f}",
        f"best_speedup_value: {best['overall_speedup']:.9f}",
        f"worst_speedup_years: {worst['years']:.3f}",
        f"worst_speedup_value: {worst['overall_speedup']:.9f}",
        f"longest_case_years: {longest['years']:.3f}",
        f"longest_case_bridge_seconds: {longest['bridge_total_seconds']:.9f}",
        f"longest_case_ias15_seconds: {longest['ias15_total_seconds']:.9f}",
        f"longest_case_speedup: {longest['overall_speedup']:.9f}",
    ]
    (outdir / "efficiency_sweep_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = load_rows(Path(args.csv))
    if not rows:
        raise SystemExit("No sweep rows found")

    plot_sweep(rows, outdir)
    write_summary(rows, outdir)
    print(f"Wrote sweep plots to {outdir}")


if __name__ == "__main__":
    main()
