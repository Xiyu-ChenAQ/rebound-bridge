from __future__ import annotations

import math
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "generalized_order"))

import plot_generalized_order as analysis


class GeneralizedOrderAnalysisTest(unittest.TestCase):
    def test_local_orders_recover_power_law(self) -> None:
        steps = [0.5, 1.0, 2.0, 4.0]
        errors = [step**4 for step in steps]

        midpoints, orders = analysis.local_orders(steps, errors)

        self.assertEqual(len(midpoints), 3)
        self.assertTrue(all(math.isclose(order, 4.0) for order in orders))

    def test_resolved_mask_rejects_reference_floor(self) -> None:
        errors = [1.0e-12, 2.0e-11, 5.0e-10]

        self.assertEqual(analysis.resolved_mask(errors, 1.0e-12), [False, True, True])


if __name__ == "__main__":
    unittest.main()
