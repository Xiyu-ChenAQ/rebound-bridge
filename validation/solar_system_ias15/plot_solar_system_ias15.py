#!/usr/bin/env python3

import argparse
import csv
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
    print(f"Wrote plots to {outdir}")


if __name__ == "__main__":
    main()
