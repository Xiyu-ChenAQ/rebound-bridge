import importlib.util
import tempfile
import unittest
from pathlib import Path


def load_plot_module():
    module_path = Path(__file__).with_name("plot_solar_system_ias15.py")
    spec = importlib.util.spec_from_file_location("plot_solar_system_ias15", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def sample_data():
    return {
        "time_yr": [0.0, 1.0, 2.0, 3.0],
        "bridge_energy_rel": [0.0, 1e-11, 2e-11, 3e-11],
        "ias15_energy_rel": [0.0, 1e-12, 2e-12, 3e-12],
        "bridge_l_rel": [0.0, 2e-11, 3e-11, 4e-11],
        "ias15_l_rel": [0.0, 2e-12, 3e-12, 4e-12],
        "bridge_em_distance_au": [0.00257, 0.00258, 0.00259, 0.00260],
        "ias15_em_distance_au": [0.00257, 0.002579, 0.002588, 0.002597],
        "bridge_laplace_deg": [0.1, 0.2, 0.25, 0.3],
        "ias15_laplace_deg": [0.1, 0.18, 0.22, 0.27],
        "energy_rel_diff": [0.0, 9e-12, 1.8e-11, 2.7e-11],
        "l_rel_diff": [0.0, 1.8e-11, 2.7e-11, 3.6e-11],
        "em_distance_diff_au": [0.0, 1e-06, 2e-06, 3e-06],
        "laplace_diff_deg": [0.0, 0.02, 0.03, 0.04],
        "earth_barycenter_error_au": [0.0, 1e-09, 2e-09, 3e-09],
        "jupiter_barycenter_error_au": [0.0, 2e-09, 3e-09, 4e-09],
    }


class PlotLongRunTest(unittest.TestCase):
    def test_long_run_plot_contract_writes_expected_files(self):
        plot_module = load_plot_module()

        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)
            plot_module.plot_long_term_stability(sample_data(), outdir)

            self.assertTrue((outdir / "bridge_longterm_stability.png").exists())
            self.assertTrue((outdir / "bridge_longterm_envelope.png").exists())


if __name__ == "__main__":
    unittest.main()
