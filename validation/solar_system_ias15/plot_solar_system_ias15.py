#!/usr/bin/env python3

import argparse
import csv
import math
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_csv(path):
    with path.open(newline="", encoding="utf-8") as fp:
        reader = csv.DictReader(fp)
        columns = {name: [] for name in reader.fieldnames or []}
        for row in reader:
            for name in columns:
                columns[name].append(float(row[name]))
    return columns


def save_fig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def abs_floor(series, floor=1e-300):
    return [max(abs(v), floor) for v in series]


def running_abs_max(series, floor=1e-300):
    out = []
    current = floor
    for value in series:
        current = max(current, abs(value))
        out.append(max(current, floor))
    return out


def running_mean(series, window):
    if not series:
        return []
    out = []
    acc = 0.0
    for i, value in enumerate(series):
        acc += value
        if i >= window:
            acc -= series[i - window]
        count = min(i + 1, window)
        out.append(acc / count)
    return out


def percentile(series, fraction):
    if not series:
        return 0.0
    values = sorted(series)
    if len(values) == 1:
        return values[0]
    position = max(0.0, min(1.0, fraction)) * (len(values) - 1)
    lo = int(math.floor(position))
    hi = int(math.ceil(position))
    if lo == hi:
        return values[lo]
    weight = position - lo
    return values[lo] * (1.0 - weight) + values[hi] * weight


def unwrap_degrees(series):
    if not series:
        return []
    out = [series[0]]
    offset = 0.0
    previous = series[0]
    for value in series[1:]:
        delta = value - previous
        if delta > 180.0:
            offset -= 360.0
        elif delta < -180.0:
            offset += 360.0
        out.append(value + offset)
        previous = value
    return out


def plot_metrics(data, outdir):
    t = data["time_yr"]

    fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)
    ax = axes[0]
    ax.semilogy(t, [abs(v) for v in data["bridge_energy_rel"]], label="Bridge")
    ax.semilogy(t, [abs(v) for v in data["ias15_energy_rel"]], label="IAS15")
    ax.set_ylabel("|dE/E0|")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    ax.semilogy(t, [abs(v) for v in data["bridge_l_rel"]], label="Bridge")
    ax.semilogy(t, [abs(v) for v in data["ias15_l_rel"]], label="IAS15")
    ax.set_ylabel("|dL/L0|")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[2]
    ax.plot(t, data["bridge_em_distance_au"], label="Bridge")
    ax.plot(t, data["ias15_em_distance_au"], label="IAS15")
    ax.set_ylabel("Earth-Moon distance [AU]")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.25)

    ax = axes[3]
    ax.plot(t, data["bridge_laplace_deg"], label="Bridge")
    ax.plot(t, data["ias15_laplace_deg"], label="IAS15")
    ax.set_ylabel("Laplace angle [deg]")
    ax.set_xlabel("Time [yr]")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.25)

    save_fig(fig, outdir / "bridge_vs_ias15_metrics.png")


def plot_residuals(data, outdir):
    t = data["time_yr"]

    fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)
    ax = axes[0]
    ax.plot(t, data["energy_rel_diff"], color="tab:red")
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_ylabel("dE/E0 diff")
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    ax.plot(t, data["l_rel_diff"], color="tab:purple")
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_ylabel("dL/L0 diff")
    ax.grid(True, alpha=0.25)

    ax = axes[2]
    ax.plot(t, data["em_distance_diff_au"], color="tab:green")
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_ylabel("Earth-Moon d diff [AU]")
    ax.grid(True, alpha=0.25)

    ax = axes[3]
    ax.plot(t, data["laplace_diff_deg"], color="tab:blue")
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_ylabel("Laplace diff [deg]")
    ax.set_xlabel("Time [yr]")
    ax.grid(True, alpha=0.25)

    save_fig(fig, outdir / "bridge_vs_ias15_residuals.png")


def plot_orbital_residuals(data, outdir):
    t = data["time_yr"]

    fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)
    ax = axes[0]
    ax.plot(t, data["earth_phase_diff_deg"], color="tab:orange")
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_ylabel("Earth phase diff [deg]")
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    ax.plot(t, data["earth_radius_diff_au"], color="tab:green")
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_ylabel("Earth radius diff [AU]")
    ax.grid(True, alpha=0.25)

    ax = axes[2]
    ax.plot(t, data["jupiter_phase_diff_deg"], color="tab:brown")
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_ylabel("Jupiter phase diff [deg]")
    ax.grid(True, alpha=0.25)

    ax = axes[3]
    ax.plot(t, data["jupiter_radius_diff_au"], color="tab:purple")
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_ylabel("Jupiter radius diff [AU]")
    ax.set_xlabel("Time [yr]")
    ax.grid(True, alpha=0.25)

    save_fig(fig, outdir / "bridge_vs_ias15_orbital_residuals.png")


def plot_long_term_stability(data, outdir):
    t = data["time_yr"]
    bridge_laplace = unwrap_degrees(data["bridge_laplace_deg"])
    bridge_laplace_shift = [value - bridge_laplace[0] for value in bridge_laplace]
    bridge_em_distance = data["bridge_em_distance_au"]
    bridge_em_distance_shift = [value - bridge_em_distance[0] for value in bridge_em_distance]

    fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)
    ax = axes[0]
    ax.semilogy(t, abs_floor(data["bridge_energy_rel"]), label="Bridge")
    ax.set_ylabel("|dE/E0|")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    ax.semilogy(t, abs_floor(data["bridge_l_rel"]), label="Bridge")
    ax.set_ylabel("|dL/L0|")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[2]
    ax.plot(t, bridge_laplace_shift, color="tab:blue")
    ax.axhline(0.0, color="0.3", lw=0.8)
    ax.set_ylabel("Laplace shift [deg]")
    ax.grid(True, alpha=0.25)

    ax = axes[3]
    ax.plot(t, bridge_em_distance, color="tab:green")
    ax.set_ylabel("Earth-Moon d [AU]")
    ax.set_xlabel("Time [yr]")
    ax.grid(True, alpha=0.25)

    save_fig(fig, outdir / "bridge_longterm_stability.png")

    fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)
    ax = axes[0]
    ax.semilogy(t, running_abs_max(data["bridge_energy_rel"]), label="Bridge")
    ax.set_ylabel("max |dE/E0|")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    ax.semilogy(t, running_abs_max(data["bridge_l_rel"]), label="Bridge")
    ax.set_ylabel("max |dL/L0|")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[2]
    ax.plot(t, running_abs_max(bridge_laplace_shift), color="tab:blue")
    ax.set_ylabel("max |Laplace shift| [deg]")
    ax.grid(True, alpha=0.25)

    ax = axes[3]
    ax.plot(t, running_abs_max(bridge_em_distance_shift), color="tab:brown")
    ax.set_ylabel("max |Earth-Moon d shift| [AU]")
    ax.set_xlabel("Time [yr]")
    ax.grid(True, alpha=0.25)

    save_fig(fig, outdir / "bridge_longterm_envelope.png")


def plot_efficiency(data, outdir):
    t = data["time_yr"]
    bridge_step_ms = [1e3 * value for value in data["bridge_step_seconds"]]
    ias15_step_ms = [1e3 * value for value in data["ias15_step_seconds"]]
    step_speedup = data["step_speedup"]
    cumulative_speedup = data["cumulative_speedup"]
    smooth_window = min(200, max(10, len(step_speedup) // 50 if step_speedup else 10))
    smooth_step_speedup = running_mean(step_speedup, smooth_window)

    fig, axes = plt.subplots(4, 1, figsize=(12, 15), sharex=True)
    ax = axes[0]
    ax.semilogy(t, abs_floor(data["bridge_step_seconds"]), label="Bridge")
    ax.semilogy(t, abs_floor(data["ias15_step_seconds"]), label="IAS15")
    ax.set_ylabel("Step time [s]")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    ax.plot(t, bridge_step_ms, label="Bridge")
    ax.plot(t, ias15_step_ms, label="IAS15")
    ax.set_ylabel("Step time [ms]")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.25)

    ax = axes[2]
    ax.plot(t, data["bridge_cumulative_seconds"], label="Bridge")
    ax.plot(t, data["ias15_cumulative_seconds"], label="IAS15")
    ax.set_ylabel("Cumulative time [s]")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.25)

    ax = axes[3]
    ax.plot(t, step_speedup, color="tab:green", alpha=0.35, label="Step speedup")
    ax.plot(t, smooth_step_speedup, color="tab:green", lw=2.0, label=f"Step speedup mean ({smooth_window})")
    ax.plot(t, cumulative_speedup, color="tab:red", lw=1.5, label="Cumulative speedup")
    ax.axhline(1.0, color="0.3", lw=0.8)
    ax.set_ylabel("IAS15 / Bridge")
    ax.set_xlabel("Time [yr]")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.25)

    save_fig(fig, outdir / "bridge_vs_ias15_efficiency.png")


def write_efficiency_summary(data, outdir):
    bridge_steps = data["bridge_step_seconds"]
    ias15_steps = data["ias15_step_seconds"]
    step_speedup = [value for value in data["step_speedup"] if math.isfinite(value) and value > 0.0]
    cumulative_speedup = [value for value in data["cumulative_speedup"] if math.isfinite(value) and value > 0.0]
    bridge_total = data["bridge_cumulative_seconds"][-1] if data["bridge_cumulative_seconds"] else 0.0
    ias15_total = data["ias15_cumulative_seconds"][-1] if data["ias15_cumulative_seconds"] else 0.0
    faster_steps = sum(1 for b, i in zip(bridge_steps, ias15_steps) if b < i)
    slower_steps = sum(1 for b, i in zip(bridge_steps, ias15_steps) if b > i)

    lines = [
        "Bridge vs IAS15 efficiency summary",
        "speedup definition: IAS15 step time / Bridge step time; values > 1 mean Bridge is faster",
        f"samples: {len(bridge_steps)}",
        f"bridge_total_seconds: {bridge_total:.9f}",
        f"ias15_total_seconds: {ias15_total:.9f}",
        f"overall_speedup: {(ias15_total / bridge_total) if bridge_total > 0.0 else 0.0:.9f}",
        f"bridge_mean_step_ms: {1e3 * statistics.fmean(bridge_steps):.6f}",
        f"ias15_mean_step_ms: {1e3 * statistics.fmean(ias15_steps):.6f}",
        f"bridge_median_step_ms: {1e3 * statistics.median(bridge_steps):.6f}",
        f"ias15_median_step_ms: {1e3 * statistics.median(ias15_steps):.6f}",
        f"bridge_p95_step_ms: {1e3 * percentile(bridge_steps, 0.95):.6f}",
        f"ias15_p95_step_ms: {1e3 * percentile(ias15_steps, 0.95):.6f}",
        f"step_speedup_median: {statistics.median(step_speedup) if step_speedup else 0.0:.9f}",
        f"step_speedup_p95: {percentile(step_speedup, 0.95) if step_speedup else 0.0:.9f}",
        f"cumulative_speedup_final: {cumulative_speedup[-1] if cumulative_speedup else 0.0:.9f}",
        f"bridge_faster_steps: {faster_steps}",
        f"ias15_faster_steps: {slower_steps}",
    ]
    (outdir / "efficiency_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    csv_path = Path(args.csv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    data = load_csv(csv_path)
    plot_metrics(data, outdir)
    plot_residuals(data, outdir)
    plot_orbital_residuals(data, outdir)
    plot_long_term_stability(data, outdir)
    plot_efficiency(data, outdir)
    write_efficiency_summary(data, outdir)
    print(f"Wrote plots to {outdir}")


if __name__ == "__main__":
    main()
