"""Build and run the numerical experiments used by the MNRAS paper."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
TARGETS = (
    "galilean_horizons_compare",
    "galilean_whfast_compare",
    "galilean_long_energy",
    "galilean_reversibility",
    "generalized_order_crossover",
    "complex_speed_ratio",
    "bridge_limiting_cases_test",
    "galilean_exomoon_survival",
)


def executable_path(build_dir: Path, config: str, target: str) -> Path:
    suffix = ".exe" if sys.platform == "win32" else ""
    configured = build_dir / config / f"{target}{suffix}"
    return configured if configured.exists() or sys.platform == "win32" else build_dir / f"{target}{suffix}"


def invoke(command: list[str | Path], dry_run: bool = False) -> None:
    rendered = [str(item) for item in command]
    print("+", subprocess.list2cmdline(rendered), flush=True)
    if not dry_run:
        subprocess.run(rendered, cwd=ROOT.parent, check=True)


def initialize_csv(path: Path, header: str, dry_run: bool) -> None:
    print("+ initialize", path, flush=True)
    if not dry_run:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(header + "\n", encoding="ascii")


def run_scaling(exe: Path, full: bool, dry_run: bool) -> None:
    out = ROOT / "results" / ("scaling" if full else "smoke") / "complex_speed_ratio.csv"
    header = (
        "n_outer,bridge_method,years,dt_outer_days,binary_ratio,whfast_dt_days,repeats,warmups,"
        "whfast_mean_seconds,bridge_mean_seconds,mean_speedup_whfast_over_bridge,"
        "whfast_best_seconds,bridge_best_seconds,best_speedup_whfast_over_bridge"
    )
    initialize_csv(out, header, dry_run)
    counts = (0, 32, 64, 128) if full else (0,)
    years = "100" if full else "0.01"
    repeats = "3" if full else "1"
    warmups = "1" if full else "0"
    for count in counts:
        invoke([exe, "--n-outer", str(count), "--years", years, "--dt-outer-days", "4", "--binary-ratio", "64", "--bridge-method", "eos-lf4-2", "--repeats", repeats, "--warmups", warmups, "--csv", out], dry_run)


def build(build_dir: Path, config: str, dry_run: bool) -> None:
    invoke(["cmake", "-S", ROOT, "-B", build_dir, "-DCMAKE_BUILD_TYPE=" + config], dry_run)
    invoke(["cmake", "--build", build_dir, "--config", config, "--target", *TARGETS], dry_run)


def run_smoke(exes: dict[str, Path], dry_run: bool) -> None:
    out = ROOT / "results" / "smoke"
    out.mkdir(parents=True, exist_ok=True)
    invoke([exes["galilean_horizons_compare"], "--method", "eos_lf4_2", "--years", "0.01", "--samples", "2", "--dt-outer-days", "0.5", "--inner-ratio", "40", "--csv", out / "galilean_compare.csv"], dry_run)
    invoke([exes["galilean_whfast_compare"], "--integrator", "whfast", "--years", "0.01", "--samples", "2", "--dt-days", "0.01", "--csv", out / "global_whfast.csv"], dry_run)
    invoke([exes["galilean_long_energy"], "--method", "bridge", "--years", "0.01", "--realizations", "1", "--outputs-per-decade", "3", "--csv", out / "long_bridge.csv"], dry_run)
    invoke([exes["galilean_reversibility"], "--csv", out / "reversibility.csv"], dry_run)
    invoke([exes["bridge_limiting_cases_test"]], dry_run)
    crossover = out / "generalized_order_crossover.csv"
    invoke([exes["generalized_order_crossover"], "--csv", crossover], dry_run)
    run_scaling(exes["complex_speed_ratio"], False, dry_run)
    for integrator in ("bridge", "ias15"):
        invoke([exes["galilean_exomoon_survival"], "--integrator", integrator, "--a-planet-au", "0.55", "--years", "0.01", "--samples", "4", "--dt-outer-days", "0.5", "--inner-ratio", "40", "--csv", out / f"exomoon_{integrator}.csv"], dry_run)


def run_speed_accuracy(exes: dict[str, Path], dry_run: bool) -> None:
    out = ROOT / "results" / "speed_accuracy"
    cost = ROOT / "results" / "cost_accuracy"
    out.mkdir(parents=True, exist_ok=True)
    cost.mkdir(parents=True, exist_ok=True)
    groups = (
        ("saba", (0.1, 0.05, 0.025, 0.02, 0.01), "galilean_saba_20yr_dt{dt}_exactcadence.csv", ()),
        ("whfast", (0.02, 0.01, 0.005, 0.002, 0.001), "galilean_whfast_20yr_dt{dt}_exactcadence.csv", ()),
        ("whfast", (0.02, 0.01, 0.005, 0.002, 0.001), "galilean_whfast_corr11c2_20yr_dt{dt}_exactcadence.csv", ("--whfast-corrector", "11", "--whfast-corrector2", "1")),
    )
    for integrator, steps, pattern, extra in groups:
        for dt in steps:
            token = f"{dt:g}".replace(".", "p")
            invoke([exes["galilean_whfast_compare"], "--integrator", integrator, "--years", "20", "--samples", "10", "--dt-days", f"{dt:g}", *extra, "--csv", out / pattern.format(dt=token)], dry_run)
    invoke([sys.executable, ROOT / "galilean" / "run_kdk_lf42_cost_accuracy.py", "--exe", exes["galilean_horizons_compare"], "--output-dir", cost, "--repeats", "5", "--years", "20", "--output-interval-days", "5"], dry_run)


def run_long_energy(exes: dict[str, Path], dry_run: bool) -> None:
    out = ROOT / "results" / "long_energy"
    out.mkdir(parents=True, exist_ok=True)
    invoke([exes["galilean_long_energy"], "--method", "bridge", "--years", "1000", "--dt-outer-days", "0.5", "--inner-ratio", "40", "--realizations", "4", "--perturbation", "0.001", "--csv", out / "calibration_lf42_bridge_1000yr_r4.csv"], dry_run)
    invoke([exes["galilean_long_energy"], "--method", "whfast", "--years", "1000", "--dt-days", "0.0085", "--whfast-corrector", "11", "--whfast-corrector2", "1", "--realizations", "4", "--perturbation", "0.001", "--csv", out / "calibration_whfast_corr11c2_dt0p0085_1000yr_r4.csv"], dry_run)
    invoke([exes["galilean_long_energy"], "--method", "bridge", "--years", "10000", "--dt-outer-days", "0.5", "--inner-ratio", "40", "--realizations", "4", "--perturbation", "0.001", "--csv", out / "ensemble_lf42_bridge_10000yr_r4.csv"], dry_run)
    invoke([exes["galilean_long_energy"], "--method", "whfast", "--years", "10000", "--dt-days", "0.0085", "--whfast-corrector", "11", "--whfast-corrector2", "1", "--realizations", "4", "--perturbation", "0.001", "--csv", out / "ensemble_whfast_corr11c2_dt0p0085_10000yr_r4.csv"], dry_run)


def exomoon_command(exe: Path, integrator: str, years: int, dt: float, phase: int, path: Path) -> list[str | Path]:
    samples = round(years * 365.25 / dt)
    return [exe, "--integrator", integrator, "--a-planet-au", "0.55", "--planet-phase-deg", str(phase), "--years", str(years), "--dt-outer-days", f"{dt:g}", "--inner-ratio", "40", "--samples", str(samples), "--csv", path]


def run_exomoon(exe: Path, dry_run: bool) -> None:
    out = ROOT / "results" / "resonance_validation"
    out.mkdir(parents=True, exist_ok=True)
    cases = [
        ("bridge", 20, 0.5, 0, "a0p55_phase0_20yr_bridge_dt0p5.csv"),
        ("bridge", 20, 0.25, 0, "a0p55_phase0_20yr_bridge_dt0p25.csv"),
        ("bridge", 20, 0.125, 0, "a0p55_phase0_20yr_bridge_dt0p125.csv"),
        ("ias15", 20, 0.5, 0, "a0p55_phase0_20yr_ias15.csv"),
        ("bridge", 20, 0.25, 120, "a0p55_phase120_20yr_bridge_dt0p25.csv"),
        ("bridge", 20, 0.25, 240, "a0p55_phase240_20yr_bridge_dt0p25.csv"),
        ("bridge", 100, 0.25, 0, "a0p55_phase0_100yr_bridge_dt0p25.csv"),
        ("ias15", 100, 0.5, 0, "a0p55_phase0_100yr_ias15.csv"),
    ]
    for integrator, years, dt, phase, filename in cases:
        invoke(exomoon_command(exe, integrator, years, dt, phase, out / filename), dry_run)


def make_plots(exes: dict[str, Path], dry_run: bool) -> None:
    scripts = [
        [sys.executable, ROOT / "galilean" / "plot_two_timestep_convergence.py", "--exe", exes["galilean_horizons_compare"], "--skip-run"],
        [sys.executable, ROOT / "generalized_order" / "plot_generalized_order.py"],
        [sys.executable, ROOT / "galilean" / "plot_20yr_speed_accuracy.py"],
        [sys.executable, ROOT / "galilean" / "plot_long_energy_ensemble.py"],
        [sys.executable, ROOT / "exomoon" / "plot_laplace_circulation_validation.py"],
    ]
    for command in scripts:
        invoke(command, dry_run)


def run_full(exes: dict[str, Path], plots_only: bool, dry_run: bool) -> None:
    generalized = ROOT / "results" / "generalized_order"
    generalized.mkdir(parents=True, exist_ok=True)
    if not plots_only:
        invoke([sys.executable, ROOT / "galilean" / "plot_two_timestep_convergence.py", "--exe", exes["galilean_horizons_compare"]], dry_run)
        invoke([exes["generalized_order_crossover"], "--csv", generalized / "generalized_order_crossover.csv"], dry_run)
        run_speed_accuracy(exes, dry_run)
        run_long_energy(exes, dry_run)
        run_scaling(exes["complex_speed_ratio"], True, dry_run)
        invoke([exes["galilean_reversibility"], "--csv", ROOT / "results" / "reversibility.csv"], dry_run)
        run_exomoon(exes["galilean_exomoon_survival"], dry_run)
    make_plots(exes, dry_run)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--build-dir", type=Path, default=ROOT / "build")
    parser.add_argument("--config", default="Release")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--plots-only", action="store_true")
    parser.add_argument("--unit-tests-only", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    build_dir = args.build_dir.resolve()
    if args.unit_tests_only:
        invoke([sys.executable, "-m", "pytest", ROOT / "tests", "-q"], args.dry_run)
        return 0
    if not args.skip_build and not args.plots_only:
        build(build_dir, args.config, args.dry_run)
    exes = {target: executable_path(build_dir, args.config, target) for target in TARGETS}
    if args.mode == "smoke":
        if args.plots_only:
            raise SystemExit("--plots-only requires --mode full")
        run_smoke(exes, args.dry_run)
        invoke([sys.executable, "-m", "pytest", ROOT / "tests", "-q"], args.dry_run)
    else:
        run_full(exes, args.plots_only, args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
