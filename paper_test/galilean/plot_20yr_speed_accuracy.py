"""Plot the 20 yr Galilean speed-accuracy comparison."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


AU_KM = 149597870.7
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "speed_accuracy"
COST_BASE = ROOT / "results" / "cost_accuracy"
STYLE = ROOT / "mnras.mplstyle"


CASES = [
    ("SABA(10,6,4)", "SABA", "#4f6f2a", "s", [
        (0.1, BASE / "galilean_saba_20yr_dt0p1_exactcadence.csv"),
        (0.05, BASE / "galilean_saba_20yr_dt0p05_exactcadence.csv"),
        (0.025, BASE / "galilean_saba_20yr_dt0p025_exactcadence.csv"),
        (0.02, BASE / "galilean_saba_20yr_dt0p02_exactcadence.csv"),
        (0.01, BASE / "galilean_saba_20yr_dt0p01_exactcadence.csv"),
    ]),
    ("naive global WHFast", "WHFast", "#8a4f18", "^", [
        (0.02, BASE / "galilean_whfast_20yr_dt0p02_exactcadence.csv"),
        (0.01, BASE / "galilean_whfast_20yr_dt0p01_exactcadence.csv"),
        (0.005, BASE / "galilean_whfast_20yr_dt0p005_exactcadence.csv"),
        (0.002, BASE / "galilean_whfast_20yr_dt0p002_exactcadence.csv"),
        (0.001, BASE / "galilean_whfast_20yr_dt0p001_exactcadence.csv"),
    ]),
    ("WHFast + correctors", "WHFast corrected", "#b77a2b", "v", [
        (0.02, BASE / "galilean_whfast_corr11c2_20yr_dt0p02_exactcadence.csv"),
        (0.01, BASE / "galilean_whfast_corr11c2_20yr_dt0p01_exactcadence.csv"),
        (0.005, BASE / "galilean_whfast_corr11c2_20yr_dt0p005_exactcadence.csv"),
        (0.002, BASE / "galilean_whfast_corr11c2_20yr_dt0p002_exactcadence.csv"),
        (0.001, BASE / "galilean_whfast_corr11c2_20yr_dt0p001_exactcadence.csv"),
    ]),
]


REPEATED_STYLES = {
    "order2": ("KDK Bridge", "#8f3b2d", "D"),
    "eos_lf4_2": ("LF(4,2) Bridge", "#245c73", "o"),
}


ANNOTATION_OFFSETS = {
    ("order2", 0.125): (8, 5),
    ("order2", 0.25): (8, 2),
    ("order2", 0.5): (8, 0),
    ("order2", 1.0): (8, 4),
    ("eos_lf4_2", 0.125): (-36, -12),
    ("eos_lf4_2", 0.25): (8, 8),
    ("eos_lf4_2", 0.5): (-34, -13),
    ("eos_lf4_2", 1.0): (8, 4),
}


def annotation_offset(method: str, dt_outer_days: float) -> tuple[int, int]:
    return ANNOTATION_OFFSETS[(method, dt_outer_days)]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"empty CSV: {path}")
    return rows


def load_repeated_summary(path: Path) -> list[dict[str, float | str]]:
    rows = read_rows(path)
    result: list[dict[str, float | str]] = []
    for row in rows:
        result.append(
            {
                "method": row["method"],
                "method_label": row["method_label"],
                "dt_outer_days": float(row["dt_outer_days"]),
                "speedup_median": float(row["speedup_median"]),
                "speedup_q25": float(row["speedup_q25"]),
                "speedup_q75": float(row["speedup_q75"]),
                "max_moon_error_km": float(row["max_moon_error_au"]) * AU_KM,
                "max_laplace_diff_deg": float(row["max_laplace_diff_deg"]),
            }
        )
    methods = {str(row["method"]) for row in result}
    if methods != set(REPEATED_STYLES):
        raise RuntimeError(f"unexpected repeated methods in {path}: {sorted(methods)}")
    return result


def summarize_case(method: str, dt_days: float, path: Path) -> dict[str, float | str]:
    rows = read_rows(path)
    last = rows[-1]
    escaped_key = "bridge_escaped" if "bridge_escaped" in last else "split_escaped"
    return {
        "method": method,
        "dt_days": dt_days,
        "samples": len(rows),
        "final_time_yr": float(last["time_yr"]),
        "speedup_vs_ias15": float(last["cumulative_speedup"]),
        "whfast_corrector": float(last.get("whfast_corrector", 0.0)),
        "whfast_corrector2": float(last.get("whfast_corrector2", 0.0)),
        "max_moon_error_au": max(float(row["max_moon_error_au"]) for row in rows),
        "max_laplace_diff_deg": max(float(row["laplace_diff_deg"]) for row in rows),
        "final_moon_error_km": float(last["max_moon_error_au"]) * AU_KM,
        "final_laplace_diff_deg": float(last["laplace_diff_deg"]),
        "escaped": float(any(row.get(escaped_key, "0") == "1" or row.get("ias15_escaped", "0") == "1" for row in rows)),
    }


def collect_summary() -> list[dict[str, float | str]]:
    summary: list[dict[str, float | str]] = []
    for method, _group, _color, _marker, cases in CASES:
        for dt_days, path in cases:
            summary.append(summarize_case(method, dt_days, path))
    return summary


def write_summary(rows: list[dict[str, float | str]], path: Path) -> None:
    fields = [
        "method",
        "dt_days",
        "samples",
        "final_time_yr",
        "speedup_vs_ias15",
        "whfast_corrector",
        "whfast_corrector2",
        "max_moon_error_au",
        "max_moon_error_km",
        "max_laplace_diff_deg",
        "final_moon_error_km",
        "final_laplace_diff_deg",
        "escaped",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["max_moon_error_km"] = float(row["max_moon_error_au"]) * AU_KM
            writer.writerow(out)


def plot(
    rows: list[dict[str, float | str]],
    repeated_rows: list[dict[str, float | str]],
    path: Path,
) -> None:
    plt.style.use(STYLE)
    fig, axes = plt.subplots(1, 2, figsize=(7.03, 3.0), constrained_layout=True)
    for method, group, color, marker, _cases in CASES:
        group_rows = sorted(
            [row for row in rows if row["method"] == method],
            key=lambda row: float(row["dt_days"]),
            reverse=True,
        )
        x = [float(row["speedup_vs_ias15"]) for row in group_rows]
        moon = [float(row["max_moon_error_au"]) * AU_KM for row in group_rows]
        laplace = [float(row["max_laplace_diff_deg"]) for row in group_rows]
        linestyle = "--" if "correctors" in method else "-"
        for ax, y in ((axes[0], moon), (axes[1], laplace)):
            ax.plot(
                x,
                y,
                marker=marker,
                color=color,
                linestyle=linestyle,
                label=group,
            )
    for method, (label, color, marker) in REPEATED_STYLES.items():
        group_rows = sorted(
            [row for row in repeated_rows if row["method"] == method],
            key=lambda row: float(row["dt_outer_days"]),
            reverse=True,
        )
        x = [float(row["speedup_median"]) for row in group_rows]
        xerr = [
            [float(row["speedup_median"]) - float(row["speedup_q25"]) for row in group_rows],
            [float(row["speedup_q75"]) - float(row["speedup_median"]) for row in group_rows],
        ]
        moon = [float(row["max_moon_error_km"]) for row in group_rows]
        laplace = [float(row["max_laplace_diff_deg"]) for row in group_rows]
        for ax, y in ((axes[0], moon), (axes[1], laplace)):
            ax.errorbar(
                x,
                y,
                xerr=xerr,
                marker=marker,
                color=color,
                linestyle="-",
                linewidth=1.6,
                markersize=4.6,
                capsize=2.0,
                elinewidth=1.1,
                capthick=1.1,
                label=label,
                zorder=4,
            )
            if ax is not axes[0]:
                continue
            endpoint_rows = (group_rows[0], group_rows[-1])
            for row, x_value, y_value in zip(group_rows, x, y):
                if row not in endpoint_rows:
                    continue
                offset = annotation_offset(method, float(row["dt_outer_days"]))
                ax.annotate(
                    f"{float(row['dt_outer_days']):g} d",
                    (x_value, y_value),
                    xytext=offset,
                    textcoords="offset points",
                    fontsize=7.0,
                    color=color,
                )
    for ax in axes:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axvline(1.0, color="0.55", linestyle=":", linewidth=1.0)
        ax.set_xlabel("speed relative to IAS15")
    axes[0].set_ylabel("position error [km]")
    axes[1].set_ylabel("Laplace-angle error [deg]")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=5)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    rows = collect_summary()
    repeated_rows = load_repeated_summary(COST_BASE / "kdk_lf42_cost_accuracy_summary.csv")
    summary = BASE / "galilean_global_baselines_20yr_exactcadence_summary.csv"
    figure = BASE / "galilean_bridge_vs_saba_whfast_20yr_exactcadence.png"
    write_summary(rows, summary)
    plot(rows, repeated_rows, figure)
    print(f"wrote {summary}")
    print(f"wrote {figure}")
    print(f"wrote {figure.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
