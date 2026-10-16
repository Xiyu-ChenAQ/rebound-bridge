"""Exercise the real paper executables' fixed-step/output schedule contract."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import subprocess
import tempfile


def run(executable: Path, output: Path, *arguments: str, succeeds: bool = True):
    command = [str(executable), *arguments, "--csv", str(output)]
    completed = subprocess.run(command, capture_output=True, text=True)
    if succeeds:
        assert completed.returncode == 0, completed.stdout + completed.stderr
        with output.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert rows
        assert all(math.isfinite(float(value)) for row in rows for value in row.values())
        return rows
    assert completed.returncode == 2, completed.stdout + completed.stderr
    assert not output.exists(), "Invalid schedules must be rejected before creating a partial CSV"
    return []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--horizons", type=Path, required=True)
    parser.add_argument("--exomoon", type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="paper-schedule-") as temporary:
        root = Path(temporary)
        # Reintroducing target - accumulated_time fails before sample 42.
        for method in ("order2", "eos_lf4_2"):
            rows = run(args.horizons, root / f"horizons-{method}.csv",
                       "--method", method, "--years", "20", "--samples", "1461",
                       "--dt-outer-days", "1", "--inner-ratio", "80")
            assert len(rows) == 1461
            assert math.isclose(float(rows[-1]["time_yr"]), 20, rel_tol=0, abs_tol=1e-12)
            assert all(abs(float(row["time_yr"]) * 365.25 - 5 * index) < 4e-10
                       for index, row in enumerate(rows, 1))
            assert all(row["bridge_ok"] == row["ias15_ok"] == "1" for row in rows)
            assert all(row["bridge_escaped"] == row["ias15_escaped"] == "0" for row in rows)

        for executable, prefix in ((args.horizons, ()),
                                   (args.exomoon, ("--a-planet-au", "0.55"))):
            for label, years, samples in (("duration", "0.011225188227241616", "1"),
                                           ("cadence", "0.010951403148528405", "3"),
                                           ("too-many", "0.010951403148528405", "9")):
                run(executable, root / f"bad-{executable.stem}-{label}.csv", *prefix,
                    "--years", years, "--samples", samples,
                    "--dt-outer-days", "0.5", "--inner-ratio", "40", succeeds=False)

        # Exomoon stores its initial state as well as each requested output.
        for integrator in ("bridge", "ias15"):
            rows = run(args.exomoon, root / f"exomoon-{integrator}.csv",
                       "--integrator", integrator, "--a-planet-au", "0.55",
                       "--years", "0.010951403148528405", "--samples", "8",
                       "--dt-outer-days", "0.5", "--inner-ratio", "40")
            assert len(rows) == 9
            assert [int(row["step"]) for row in rows] == list(range(9))
            assert float(rows[0]["time_yr"]) == 0
            assert abs(float(rows[-1]["time_yr"]) * 365.25 - 4) < 1e-12
    print("Paper harness schedule regressions passed")


if __name__ == "__main__":
    main()
