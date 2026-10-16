"""Reproduce the seven manuscript figures, or run short executable checks.

--dry-run prints a complete command plan without creating files.
Full experiments are expensive; --experiments selects individual figures.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parents[1]
TARGETS = (
    "complex_speed_ratio", "galilean_horizons_compare", "galilean_whfast_compare",
    "galilean_long_energy", "galilean_reversibility", "generalized_order_crossover",
    "bridge_limiting_cases_test", "galilean_exomoon_survival", "multiple_subsystems",
)
FIGURES = {
    "convergence": "galilean_two_timestep_convergence.pdf",
    "generalized_order": "generalized_order_crossover_column.pdf",
    "speed_accuracy": "galilean_bridge_vs_saba_whfast_20yr_exactcadence_column.pdf",
    "scaling": "runtime_scaling_column.pdf",
    "multiple_subsystems": "multiple_subsystems.pdf",
    "long_energy": "galilean_long_ensemble_10000yr_r4_column.pdf",
    "exomoon": "laplace_circulation_validation_column.pdf",
}
EXPERIMENTS = (*FIGURES, "reversibility")


def executable_path(build_dir: Path, config: str, target: str) -> Path:
    suffix = ".exe" if sys.platform == "win32" else ""
    configured = build_dir / config / f"{target}{suffix}"
    single = build_dir / f"{target}{suffix}"
    if configured.exists():
        return configured
    if single.exists():
        return single
    return configured if sys.platform == "win32" else single


def exomoon_command(exe: Path, integrator: str, years: float, dt: float,
                    phase: int, path: Path) -> list[str | Path]:
    samples = round(years * 365.25 / dt)
    return [exe, "--integrator", integrator, "--a-planet-au", "0.55",
            "--planet-phase-deg", str(phase), "--years", f"{years:.17g}",
            "--dt-outer-days", f"{dt:g}", "--inner-ratio", "40",
            "--samples", str(samples), "--csv", path]


def resolve_executables(argv: list[str], build_dir: Path, config: str) -> list[str]:
    # The generator may choose a different layout from the pre-build plan.
    suffix = ".exe" if sys.platform == "win32" else ""
    replacements = {}
    for name in TARGETS:
        actual = str(executable_path(build_dir, config, name))
        replacements[str(build_dir / f"{name}{suffix}")] = actual
        replacements[str(build_dir / config / f"{name}{suffix}")] = actual
    return [replacements.get(arg, arg) for arg in argv]


def make_plan(args: argparse.Namespace) -> dict:
    out, build_dir = args.output_dir.resolve(), args.build_dir.resolve()
    exes = {name: executable_path(build_dir, args.config, name) for name in TARGETS}
    commands: list[dict] = []
    selected = args.experiments or list(EXPERIMENTS)
    full = args.mode == "full"
    figures = out / "figures"

    def add(experiment: str, kind: str, *command: str | Path) -> None:
        commands.append({"experiment": experiment, "kind": kind,
                         "argv": [str(value) for value in command]})

    def script(experiment: str, kind: str, path: str, *options: str | Path) -> None:
        add(experiment, kind, sys.executable, ROOT / path, *options)

    if not args.skip_build and not args.plots_only:
        add("build", "build", "cmake", "-S", ROOT, "-B", build_dir, "-DCMAKE_BUILD_TYPE=" + args.config)
        add("build", "build", "cmake", "--build", build_dir, "--config", args.config, "--target", *TARGETS)
    if not args.plots_only:
        add("checks", "test", "ctest", "--test-dir", build_dir, "-C", args.config, "--output-on-failure")

    for experiment in selected:
        if experiment == "convergence":
            data = out / "two_timestep_convergence"
            if not args.plots_only:
                if full:
                    script(experiment, "integration", "galilean/plot_two_timestep_convergence.py",
                           "--exe", exes["galilean_horizons_compare"], "--output-dir", data)
                else:
                    for method in ("order2", "eos_lf4_2"):
                        add(experiment, "integration", exes["galilean_horizons_compare"],
                            "--method", method, "--years", f"{8 / 365.25:.17g}",
                            "--samples", "8", "--dt-outer-days", "0.5", "--inner-ratio", "40",
                            "--csv", data / f"smoke_{method}.csv")
            if full:
                script(experiment, "plot", "plot_convergence_figures.py", "--only", "galilean",
                       "--galilean-summary", data / "galilean_two_timestep_convergence_summary.csv",
                       "--output-dir", figures)
        elif experiment == "generalized_order":
            csv = out / experiment / "generalized_order_crossover.csv"
            if not args.plots_only:
                add(experiment, "integration", exes["generalized_order_crossover"], "--csv", csv)
            if full:
                script(experiment, "plot", "plot_convergence_figures.py", "--only", experiment,
                       "--crossover-csv", csv, "--output-dir", figures)
        elif experiment == "speed_accuracy":
            data, cost = out / experiment, out / "cost_accuracy"
            if not args.plots_only:
                groups = (
                    ("saba", (0.1, 0.05, 0.025, 0.02, 0.01), "saba", ()),
                    ("whfast", (0.02, 0.01, 0.005, 0.002, 0.001), "whfast", ()),
                    ("whfast", (0.02, 0.01, 0.005, 0.002, 0.001), "whfast_corr11c2",
                     ("--whfast-corrector", "11", "--whfast-corrector2", "1")),
                )
                for integrator, steps, tag, extra in groups:
                    for dt in steps if full else steps[:1]:
                        token = f"{dt:g}".replace(".", "p")
                        name = f"galilean_{tag}_20yr_dt{token}_exactcadence.csv" if full else f"smoke_{tag}.csv"
                        add(experiment, "integration", exes["galilean_whfast_compare"],
                            "--integrator", integrator, "--years", "20" if full else f"{8/365.25:.17g}",
                            "--samples", "1461" if full else "8", "--dt-days", str(dt), *extra, "--csv", data / name)
                if full:
                    script(experiment, "integration", "galilean/run_kdk_lf42_cost_accuracy.py",
                           "--exe", exes["galilean_horizons_compare"], "--output-dir", cost,
                           "--repeats", "5", "--years", "20", "--output-interval-days", "5",
                           "--build-config", args.config)
            if full:
                script(experiment, "plot", "galilean/plot_20yr_speed_accuracy.py",
                       "--data-dir", data, "--cost-dir", cost, "--output-dir", figures)
        elif experiment == "scaling":
            csv = out / experiment / "complex_speed_ratio.csv"
            if not args.plots_only:
                for count in (0, 4, 8, 16, 32, 64, 128) if full else (0, 4):
                    add(experiment, "integration", exes["complex_speed_ratio"], "--n-outer", str(count),
                        "--years", "100" if full else f"{8/365.25:.17g}",
                        "--dt-outer-days", "4", "--binary-ratio", "64", "--bridge-method", "eos-lf4-2",
                        "--repeats", "3" if full else "1", "--warmups", "1" if full else "0", "--csv", csv)
            if full:
                script(experiment, "plot", "scaling/plot_runtime_scaling.py",
                       "--csv", csv, "--output", figures / "runtime_scaling_column")
        elif experiment == "multiple_subsystems":
            data = out / experiment
            if not args.plots_only:
                script(experiment, "integration", "multiple_subsystems/run_experiment.py",
                       "--exe", exes["multiple_subsystems"], "--output", data, *([] if full else ["--pilot"]))
            if full:
                script(experiment, "plot", "plot_convergence_figures.py", "--only", experiment,
                       "--multiple-data", data, "--output-dir", figures)
        elif experiment == "long_energy":
            data = out / experiment
            if not args.plots_only:
                for method in ("bridge", "whfast"):
                    tag = "lf42_bridge" if method == "bridge" else "whfast_corr11c2_dt0p0085"
                    name = f"ensemble_{tag}_10000yr_r4.csv" if full else f"smoke_{method}.csv"
                    add(experiment, "integration", exes["galilean_long_energy"], "--method", method,
                        "--years", "10000" if full else f"{8/365.25:.17g}",
                        *(["--dt-outer-days", "0.5", "--inner-ratio", "40"] if method == "bridge"
                          else ["--dt-days", "0.0085", "--whfast-corrector", "11", "--whfast-corrector2", "1"]),
                        "--realizations", "4" if full else "1", "--perturbation", "0.001", "--csv", data / name)
            if full:
                script(experiment, "plot", "galilean/plot_long_energy_ensemble.py",
                       "--data-dir", data, "--output-dir", figures)
        elif experiment == "exomoon":
            data = out / "resonance_validation"
            cases = [
                ("bridge", 20, 0.5, 0), ("bridge", 20, 0.25, 0), ("bridge", 20, 0.125, 0),
                ("ias15", 20, 0.5, 0), ("bridge", 20, 0.25, 120), ("bridge", 20, 0.25, 240),
                ("bridge", 100, 0.25, 0), ("ias15", 100, 0.5, 0),
            ] if full else [(method, 8/365.25, 0.5, 0) for method in ("bridge", "ias15")]
            if not args.plots_only:
                for integrator, years, dt, phase in cases:
                    tag = f"{dt:g}".replace(".", "p")
                    name = (f"a0p55_phase{phase}_{years}yr_{integrator}"
                            + (f"_dt{tag}" if integrator == "bridge" else "") + ".csv") if full else f"smoke_{integrator}.csv"
                    add(experiment, "integration", *exomoon_command(
                        exes["galilean_exomoon_survival"], integrator, years, dt, phase, data / name))
            if full:
                script(experiment, "plot", "exomoon/plot_laplace_circulation_validation.py",
                       "--data-dir", data, "--output-dir", figures)
        elif experiment == "reversibility" and not args.plots_only:
            add(experiment, "integration", exes["galilean_reversibility"], "--csv", out / "reversibility.csv")
    return {"mode": args.mode, "experiments": selected, "output_dir": str(out),
            "build_dir": str(build_dir), "config": args.config, "commands": commands,
            "figures": [str(figures / FIGURES[name]) for name in selected if name in FIGURES] if full else [],
            "scope": "Manuscript configurations" if full else "Short execution checks; not manuscript results"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_hashes() -> dict[str, str]:
    files = [REPOSITORY_ROOT / name for name in ("CMakeLists.txt", "pyproject.toml", "MANIFEST.in")]
    for directory in ("src", "rebound/src", "tests/paper", "cmake"):
        files.extend(path for path in (REPOSITORY_ROOT / directory).rglob("*")
                     if path.is_file() and path.suffix in {".c", ".h", ".py", ".md", ".txt", ".cmake", ".in", ".mplstyle"}
                     and "results" not in path.relative_to(REPOSITORY_ROOT).parts)
    return {path.relative_to(REPOSITORY_ROOT).as_posix(): sha256(path)
            for path in sorted(set(files)) if path.is_file()}


def execute(plan: dict, plots_only: bool) -> int:
    out = Path(plan["output_dir"])
    if not plots_only and out.exists() and any(out.iterdir()):
        raise SystemExit(f"Output directory is not empty: {out}. Choose a fresh --output-dir; use --plots-only to redraw existing data.")
    if plots_only and not out.is_dir():
        raise SystemExit(f"Existing results required for --plots-only: {out}")
    prior_manifest = out / "run_manifest.json"
    if plots_only and prior_manifest.is_file():
        prior = json.loads(prior_manifest.read_text(encoding="utf-8"))
        if prior.get("mode") == "smoke":
            raise SystemExit("Cannot redraw manuscript figures from smoke results; use full experiment data.")
    out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    logs = out / "logs" / stamp
    logs.mkdir(parents=True)
    manifest_path = out / (f"plot_manifest_{stamp}.json" if plots_only else "run_manifest.json")
    manifest = {**plan, "started_utc": datetime.now(timezone.utc).isoformat(),
                "python": sys.version, "platform": platform.platform(),
                "source_sha256": source_hashes(), "status": "running", "executions": []}

    def save() -> None:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    save()
    try:
        for index, command in enumerate(plan["commands"], 1):
            argv = resolve_executables(command["argv"], Path(plan["build_dir"]), plan["config"])
            if "--csv" in argv:
                Path(argv[argv.index("--csv") + 1]).parent.mkdir(parents=True, exist_ok=True)
            log_path = logs / f"{index:03d}-{command['experiment']}.log"
            record = {**command, "argv": argv, "log": str(log_path),
                      "started_utc": datetime.now(timezone.utc).isoformat()}
            manifest["executions"].append(record)
            save()
            print(f"[{index}/{len(plan['commands'])}] {command['experiment']}: {subprocess.list2cmdline(argv)}", flush=True)
            started = time.monotonic()
            with log_path.open("w", encoding="utf-8") as log:
                process = subprocess.run(argv, cwd=REPOSITORY_ROOT, stdout=log, stderr=subprocess.STDOUT)
            record.update(returncode=process.returncode, elapsed_seconds=time.monotonic()-started)
            save()
            if process.returncode:
                raise RuntimeError(f"Command failed with exit code {process.returncode}; see {log_path}")
        missing = [path for path in plan["figures"] if not Path(path).is_file()]
        if missing:
            raise RuntimeError("Expected figures were not written: " + ", ".join(missing))
        manifest["status"] = "complete"
    except BaseException as exc:
        manifest["status"] = "failed"
        manifest["error"] = str(exc)
        raise
    finally:
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["output_sha256"] = {path.relative_to(out).as_posix(): sha256(path)
                                     for path in sorted(out.rglob("*")) if path.is_file() and path != manifest_path}
        manifest["executable_sha256"] = {name: sha256(path) for name in TARGETS
            if (path := executable_path(Path(plan["build_dir"]), plan["config"], name)).is_file()}
        save()
    print(f"Completed {plan['mode']} run. Record: {manifest_path}", flush=True)
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--build-dir", type=Path, default=REPOSITORY_ROOT / "build" / "paper")
    parser.add_argument("--output-dir", type=Path, help="fresh results directory; defaults to results/<mode>")
    parser.add_argument("--config", default="Release")
    parser.add_argument("--experiments", nargs="+", choices=EXPERIMENTS)
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--plots-only", action="store_true")
    parser.add_argument("--unit-tests-only", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    if args.plots_only and args.mode != "full":
        parser.error("--plots-only requires --mode full")
    if args.output_dir is None:
        args.output_dir = ROOT / "results" / args.mode
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.unit_tests_only:
        command = [sys.executable, "-m", "pytest", str(ROOT / "checks"), "-q"]
        if args.dry_run:
            print(json.dumps({"commands": [command]}, indent=2))
            return 0
        return subprocess.run(command, cwd=REPOSITORY_ROOT).returncode
    plan = make_plan(args)
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return 0
    return execute(plan, args.plots_only)


if __name__ == "__main__":
    raise SystemExit(main())
