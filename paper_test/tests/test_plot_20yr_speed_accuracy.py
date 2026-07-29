from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "galilean"))

from plot_20yr_speed_accuracy import (
    CASES,
    REPEATED_STYLES,
    annotation_offset,
    load_repeated_summary,
    plot,
)


class RepeatedSummaryPlotTests(unittest.TestCase):
    def test_plot_writes_raster_preview_and_vector_paper_figure(self) -> None:
        baseline_rows = [
            {
                "method": method,
                "dt_days": 0.1,
                "speedup_vs_ias15": 0.5 + index,
                "max_moon_error_au": 1.0e-6 * (index + 1),
                "max_laplace_diff_deg": 1.0e-4 * (index + 1),
            }
            for index, (method, _group, _color, _marker, _cases) in enumerate(CASES)
        ]
        repeated_rows = [
            {
                "method": method,
                "dt_outer_days": 0.125,
                "speedup_median": 2.0 + index,
                "speedup_q25": 1.9 + index,
                "speedup_q75": 2.1 + index,
                "max_moon_error_km": 100.0 * (index + 1),
                "max_laplace_diff_deg": 1.0e-3 * (index + 1),
            }
            for index, method in enumerate(REPEATED_STYLES)
        ]

        with tempfile.TemporaryDirectory() as directory:
            png = Path(directory) / "speed_accuracy.png"
            plot(baseline_rows, repeated_rows, png)

            self.assertTrue(png.exists())
            self.assertTrue(png.with_suffix(".pdf").exists())

    def test_annotation_offsets_separate_the_lf_error_floor(self) -> None:
        self.assertNotEqual(
            annotation_offset("eos_lf4_2", 0.125),
            annotation_offset("eos_lf4_2", 0.25),
        )
        self.assertNotEqual(
            annotation_offset("eos_lf4_2", 0.5),
            annotation_offset("order2", 0.125),
        )

    def test_loader_preserves_timing_intervals_and_converts_distance(self) -> None:
        fields = [
            "method",
            "method_label",
            "dt_outer_days",
            "speedup_median",
            "speedup_q25",
            "speedup_q75",
            "max_moon_error_au",
            "max_laplace_diff_deg",
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "summary.csv"
            with path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                writer.writerow(
                    {
                        "method": "order2",
                        "method_label": "KDK / order2",
                        "dt_outer_days": 0.5,
                        "speedup_median": 2.0,
                        "speedup_q25": 1.8,
                        "speedup_q75": 2.2,
                        "max_moon_error_au": 1.0e-6,
                        "max_laplace_diff_deg": 3.0e-4,
                    }
                )
                writer.writerow(
                    {
                        "method": "eos_lf4_2",
                        "method_label": "EOS LF4_2",
                        "dt_outer_days": 0.5,
                        "speedup_median": 1.9,
                        "speedup_q25": 1.7,
                        "speedup_q75": 2.0,
                        "max_moon_error_au": 2.0e-7,
                        "max_laplace_diff_deg": 5.0e-5,
                    }
                )

            rows = load_repeated_summary(path)

        self.assertEqual({row["method"] for row in rows}, {"order2", "eos_lf4_2"})
        self.assertEqual(rows[0]["speedup_q25"], 1.8)
        self.assertEqual(rows[0]["speedup_q75"], 2.2)
        self.assertAlmostEqual(rows[0]["max_moon_error_km"], 149.5978707)
        self.assertEqual(rows[1]["max_laplace_diff_deg"], 5.0e-5)


if __name__ == "__main__":
    unittest.main()
