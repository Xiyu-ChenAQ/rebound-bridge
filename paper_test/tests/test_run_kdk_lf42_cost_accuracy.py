from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "galilean"))

from run_kdk_lf42_cost_accuracy import (
    aggregate_case,
    build_cases,
    build_run_schedule,
    compute_work_counts,
    inclusive_quantile,
    validate_configuration,
)


class ScheduleAndWorkCountTests(unittest.TestCase):
    def test_build_cases_keeps_the_internal_timestep_fixed(self) -> None:
        cases = build_cases()

        self.assertEqual(len(cases), 8)
        self.assertEqual(
            [(case.dt_outer_days, case.inner_ratio) for case in cases[:4]],
            [(0.125, 10), (0.25, 20), (0.5, 40), (1.0, 80)],
        )
        self.assertTrue(all(case.dt_inner_days == 0.0125 for case in cases))

    def test_kdk_work_counts_include_output_call_boundaries(self) -> None:
        counts = compute_work_counts("order2", 7305.0, 1461, 0.5, 0.0125)

        self.assertEqual(counts.outer_steps, 14610)
        self.assertEqual(counts.flow_segments, 14610)
        self.assertEqual(counts.cross_kicks, 16071)
        self.assertEqual(counts.internal_steps, 584400)

    def test_lf42_work_counts_include_equal_segment_substeps(self) -> None:
        counts = compute_work_counts("eos_lf4_2", 7305.0, 1461, 0.5, 0.0125)

        self.assertEqual(counts.outer_steps, 14610)
        self.assertEqual(counts.flow_segments, 43830)
        self.assertEqual(counts.cross_kicks, 29220)
        self.assertEqual(counts.internal_steps, 613620)

    def test_inclusive_quantile_interpolates_at_n_minus_one(self) -> None:
        values = [1.0, 2.0, 3.0, 4.0, 5.0]

        self.assertEqual(inclusive_quantile(values, 0.25), 2.0)
        self.assertEqual(inclusive_quantile(values, 0.75), 4.0)

    def test_run_schedule_is_complete_and_deterministic_per_repeat(self) -> None:
        first = build_run_schedule(repeats=2, seed=1234)
        second = build_run_schedule(repeats=2, seed=1234)

        self.assertEqual(first, second)
        self.assertEqual(len(first), 16)
        for repeat in (1, 2):
            round_cases = [item.case for item in first if item.repeat == repeat]
            self.assertCountEqual(round_cases, build_cases())
            self.assertEqual(
                [item.order for item in first if item.repeat == repeat],
                list(range(1, 9)),
            )

    def test_configuration_requires_exact_common_output_targets(self) -> None:
        self.assertEqual(validate_configuration(20.0, 5.0, build_cases()), 1461)

        with self.assertRaisesRegex(ValueError, "output intervals"):
            validate_configuration(20.0, 7.0, build_cases())


class AggregationTests(unittest.TestCase):
    def _write_run(
        self,
        path: Path,
        method_seconds: float,
        ias15_seconds: float,
        moon_error: float = 2.0e-6,
        laplace_error: float = 3.0e-4,
    ) -> None:
        fields = [
            "time_yr",
            "bridge_cumulative_seconds",
            "ias15_cumulative_seconds",
            "cumulative_speedup",
            "bridge_escaped",
            "ias15_escaped",
            "max_moon_error_au",
            "laplace_diff_deg",
        ]
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerow(
                {
                    "time_yr": 10.0,
                    "bridge_cumulative_seconds": method_seconds / 2.0,
                    "ias15_cumulative_seconds": ias15_seconds / 2.0,
                    "cumulative_speedup": ias15_seconds / method_seconds,
                    "bridge_escaped": 0,
                    "ias15_escaped": 0,
                    "max_moon_error_au": moon_error / 2.0,
                    "laplace_diff_deg": laplace_error / 2.0,
                }
            )
            writer.writerow(
                {
                    "time_yr": 20.0,
                    "bridge_cumulative_seconds": method_seconds,
                    "ias15_cumulative_seconds": ias15_seconds,
                    "cumulative_speedup": ias15_seconds / method_seconds,
                    "bridge_escaped": 0,
                    "ias15_escaped": 0,
                    "max_moon_error_au": moon_error,
                    "laplace_diff_deg": laplace_error,
                }
            )

    def test_aggregate_case_reports_median_and_inclusive_quartiles(self) -> None:
        case = build_cases()[0]
        with tempfile.TemporaryDirectory() as directory:
            paths: list[Path] = []
            for repeat, method_seconds in enumerate((4.0, 5.0, 6.0, 7.0, 8.0)):
                path = Path(directory) / f"run_{repeat}.csv"
                self._write_run(path, method_seconds, 12.0)
                paths.append(path)

            summary = aggregate_case(case, paths, expected_years=20.0, expected_samples=2)

        self.assertEqual(summary.repeats, 5)
        self.assertEqual(summary.method_seconds_median, 6.0)
        self.assertEqual(summary.method_seconds_q25, 5.0)
        self.assertEqual(summary.method_seconds_q75, 7.0)
        self.assertEqual(summary.ias15_seconds_median, 12.0)
        self.assertEqual(summary.speedup_median, 2.0)
        self.assertAlmostEqual(summary.speedup_q25, 12.0 / 7.0)
        self.assertAlmostEqual(summary.speedup_q75, 12.0 / 5.0)
        self.assertEqual(summary.max_moon_error_au, 2.0e-6)
        self.assertEqual(summary.max_laplace_diff_deg, 3.0e-4)

    def test_aggregate_case_rejects_repeat_dependent_accuracy(self) -> None:
        case = build_cases()[0]
        with tempfile.TemporaryDirectory() as directory:
            paths: list[Path] = []
            for repeat in range(5):
                path = Path(directory) / f"run_{repeat}.csv"
                error = 2.0e-6 if repeat < 4 else 2.1e-6
                self._write_run(path, 5.0 + repeat, 12.0, moon_error=error)
                paths.append(path)

            with self.assertRaisesRegex(ValueError, "repeat-dependent"):
                aggregate_case(case, paths, expected_years=20.0, expected_samples=2)


if __name__ == "__main__":
    unittest.main()
