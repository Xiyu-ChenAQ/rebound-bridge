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


def plot_host_errors(data, outdir):
    t = data["time_yr"]

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    ax = axes[0]
    ax.semilogy(t, data["earth_barycenter_error_au"], color="tab:orange")
    ax.set_ylabel("Earth-Moon host error [AU]")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    ax.semilogy(t, data["jupiter_barycenter_error_au"], color="tab:brown")
    ax.set_ylabel("Jupiter host error [AU]")
    ax.set_xlabel("Time [yr]")
    ax.grid(True, which="both", alpha=0.25)

    save_fig(fig, outdir / "bridge_host_errors.png")


def plot_long_term_stability(data, outdir):
    t = data["time_yr"]

    fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)
    ax = axes[0]
    ax.semilogy(t, abs_floor(data["bridge_energy_rel"]), label="Bridge")
    ax.semilogy(t, abs_floor(data["ias15_energy_rel"]), label="IAS15")
    ax.set_ylabel("|dE/E0|")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    ax.semilogy(t, abs_floor(data["bridge_l_rel"]), label="Bridge")
    ax.semilogy(t, abs_floor(data["ias15_l_rel"]), label="IAS15")
    ax.set_ylabel("|dL/L0|")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[2]
    ax.semilogy(t, abs_floor(data["laplace_diff_deg"]), color="tab:blue")
    ax.set_ylabel("|Laplace diff| [deg]")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[3]
    ax.semilogy(t, abs_floor(data["earth_barycenter_error_au"]), label="Earth-Moon host")
    ax.semilogy(t, abs_floor(data["jupiter_barycenter_error_au"]), label="Jupiter host")
    ax.set_ylabel("Host error [AU]")
    ax.set_xlabel("Time [yr]")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)

    save_fig(fig, outdir / "bridge_longterm_stability.png")

    fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)
    ax = axes[0]
    ax.semilogy(t, running_abs_max(data["bridge_energy_rel"]), label="Bridge")
    ax.semilogy(t, running_abs_max(data["ias15_energy_rel"]), label="IAS15")
    ax.set_ylabel("max |dE/E0|")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[1]
    ax.semilogy(t, running_abs_max(data["bridge_l_rel"]), label="Bridge")
    ax.semilogy(t, running_abs_max(data["ias15_l_rel"]), label="IAS15")
    ax.set_ylabel("max |dL/L0|")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[2]
    ax.semilogy(t, running_abs_max(data["laplace_diff_deg"]), color="tab:blue")
    ax.set_ylabel("max |Laplace diff| [deg]")
    ax.grid(True, which="both", alpha=0.25)

    ax = axes[3]
    host_error = [max(a, b) for a, b in zip(
        abs_floor(data["earth_barycenter_error_au"]),
        abs_floor(data["jupiter_barycenter_error_au"]),
    )]
    ax.semilogy(t, running_abs_max(host_error), color="tab:brown")
    ax.set_ylabel("max host error [AU]")
    ax.set_xlabel("Time [yr]")
    ax.grid(True, which="both", alpha=0.25)

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
    plot_host_errors(data, outdir)
    plot_long_term_stability(data, outdir)
    print(f"Wrote plots to {outdir}")


if __name__ == "__main__":
    main()
