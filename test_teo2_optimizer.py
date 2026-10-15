"""Numerical/physical regression checks, not experimental validation."""
import math
import json
import subprocess
import sys
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np

from teo2_optimizer import (Material, Model, Settings, acoustic_modes,
                            christoffel, impermeability, indices, response,
                            rotational_impermeability, flux_projection, two_wave_efficiency,
                            optimize, optimize_design)
from fit_dispersion import fit, WAVELENGTH_UM, OBSERVED, UNCERTAINTY


class PhysicsChecks(unittest.TestCase):
    def setUp(self):
        self.material = Material()
        self.settings = Settings(632.8, rf_min_mhz=0.01, rf_max_mhz=1000,
                                 rf_budget_w=0.01, optical_transmission=0.95)
        self.model = Model(self.settings)
        self.record = self.model.record(0, 15, 0, "o", "e", 1)
        self.assertIsNotNone(self.record)
        self.report = dict(inputs=asdict(self.settings), material=asdict(self.material), optimum=self.record)

    def test_material_stability_and_published_velocity_scale(self):
        self.material.validate()
        q = np.array([1, 1, 0])/math.sqrt(2)
        velocities, u = acoustic_modes(q, self.material)
        expected = math.sqrt((self.material.c11-self.material.c12)*1e9/(2*self.material.rho))
        self.assertAlmostEqual(velocities[0], expected, places=8)
        self.assertAlmostEqual(velocities[0], 617, delta=4)  # literature benchmark, not exact same data set
        for j in range(3):
            np.testing.assert_allclose(christoffel(q, self.material)@u[:, j],
                                       self.material.rho*velocities[j]**2*u[:, j], rtol=1e-12, atol=1e-4)

    def test_shear_convention_and_slow_shear_coupling(self):
        b = impermeability(np.array([1., 0, 0]), np.array([0., 1, 0]), self.material)
        self.assertAlmostEqual(b[0, 1], self.material.p66)
        q, u = np.array([1., 1, 0])/math.sqrt(2), np.array([1., -1, 0])/math.sqrt(2)
        b = impermeability(q, u, self.material)
        self.assertAlmostEqual(b[0, 0], (self.material.p11-self.material.p12)/2)
        self.assertAlmostEqual(b[1, 1], -(self.material.p11-self.material.p12)/2)

    def test_dispersion_and_extraordinary_surface_limits(self):
        no, ne = indices(632.8)
        self.assertAlmostEqual(no, 2.2597, delta=0.001)
        self.assertAlmostEqual(ne, 2.4119, delta=0.001)
        n, _ = self.model.optical(self.model.direction([0, 90]), "e")
        np.testing.assert_allclose(n, [no, ne], atol=1e-12)
        with self.assertRaises(ValueError):
            indices(1064)
        self.assertTrue(all(np.isfinite(indices(1064, True))))

    def test_primary_index_data_and_reproducible_fit(self):
        fits = fit()
        for branch in fits["branches"].values():
            self.assertLess(branch["max_abs_index_residual"], 0.0005)
        for i, w in enumerate(WAVELENGTH_UM):
            values = indices(w*1000)
            for j, pol in enumerate("oe"):
                self.assertLess(abs(values[j]-OBSERVED[pol][i]), UNCERTAINTY[pol][i])

    def test_phase_matching_and_literature_figure_merit_scale(self):
        r = self.model.branch(0, [15], 0, "o", "e", 1)
        self.assertLess(r["residual"][0], 1e-12)
        self.assertGreater(r["m2"][0]/1e-15, 650)
        self.assertLess(r["m2"][0]/1e-15, 850)
        q = r["q"]
        bi = r["ni"][0]*self.model.direction([15])[0]
        bd = bi + r["Q"][0]*q
        np.testing.assert_allclose(bd-bi, r["Q"][0]*q, atol=1e-14)

    def test_response_matches_on_resonance_analytic_result(self):
        eta, mismatch = response(self.report)
        self.assertAlmostEqual(float(eta), self.record["throughput_efficiency"], places=10)
        self.assertLess(abs(float(mismatch)), 1e-6)

    def test_power_zero_first_peak_overdrive_and_small_signal(self):
        peak = self.record["rf_first_peak_w"]
        powers = np.array([0, 1e-8, peak, 4*peak])
        eta, _ = response(self.report, rf_power_w=powers)
        self.assertEqual(eta[0], 0)
        self.assertAlmostEqual(eta[1]/1e-8, self.record["gain_per_rf_w"]*0.95, delta=1e-4)
        self.assertAlmostEqual(eta[2], 0.95, places=10)
        self.assertLess(eta[3], 1e-20)

    def test_detuning_and_masked_near_axis_region(self):
        f = self.record["rf_mhz"]
        eta, mismatch = response(self.report, rf_mhz=np.array([f, 1.5*f]), rf_power_w=0.01)
        self.assertLess(abs(mismatch[0]), 1e-6)
        self.assertGreater(abs(mismatch[1]), 100)
        self.assertLess(eta[1], eta[0])
        masked, _ = response(self.report, alignment_offset_deg=-13)
        self.assertTrue(np.isnan(masked))

    def test_bad_inputs_rejected(self):
        with self.assertRaises(ValueError):
            Model(Settings(632.8, conversion=0))
        with self.assertRaises(ValueError):
            Model(Settings(float("nan")))
        with self.assertRaises(ValueError):
            response(self.report, rf_power_w=-1)

    def test_rotation_matches_finite_rotated_tensor_and_isotropic_limit(self):
        from scipy.linalg import expm
        q, u = np.array([0., 0, 1]), np.array([1., 0, 0])
        no, ne = indices(632.8)
        h = np.outer(u, q)
        w = (h-h.T)/2
        b = np.diag([1/no**2, 1/no**2, 1/ne**2])
        step = 1e-6
        rp, rm = expm(step*w), expm(-step*w)
        numerical = (rp@b@rp.T-rm@b@rm.T)/(2*step)
        np.testing.assert_allclose(numerical, rotational_impermeability(q, u, no, ne), atol=1e-10)
        np.testing.assert_allclose(rotational_impermeability(q, u, no, no), np.zeros((3, 3)), atol=1e-15)

    def test_poynting_flux_isotropic_oblique_limit(self):
        n = 2.3
        theta = math.radians(25)
        s = np.array([math.sin(theta), 0, math.cos(theta)])
        d = np.array([0., 1, 0])
        j = flux_projection(s, d, n, n, np.array([0., 0, 1]))[0]
        self.assertAlmostEqual(j, math.cos(theta)/n**3, places=12)

    def test_chang_two_wave_closed_form(self):
        # Compare to a independently evolved 2x2 coupled-amplitude system.
        from scipy.linalg import expm
        length, kappa, delta = 0.005, 180.0, 300.0
        hamiltonian = np.array([[delta/2, kappa], [kappa, -delta/2]])
        amplitudes = expm(-1j*hamiltonian*length)@np.array([1., 0.])
        self.assertAlmostEqual(float(two_wave_efficiency(kappa**2, delta, length)),
                               abs(amplitudes[1])**2, places=12)
        self.assertAlmostEqual(sum(abs(amplitudes)**2), 1, places=12)

    def test_parallel_tangents_and_first_order_angular_stationarity(self):
        settings = Settings(632.8, design_mode="wide-aperture", rf_min_mhz=0.001,
                            rf_max_mhz=1000, rf_budget_w=0.01)
        model = Model(settings)
        a, th, pd, root, r = model.wide_geometry(30, 0, "o")
        self.assertLess(abs(r["parallel_cross"][0]), 1e-12)
        rec = model.record(a, th, 0, "o", pd, root)
        self.assertIsNotNone(rec)
        report = dict(inputs=asdict(settings), material=asdict(self.material), optimum=rec)
        _, mismatch = response(report, alignment_offset_deg=np.array([-0.001, 0, 0.001]))
        slope = (mismatch[2]-mismatch[0])/0.002
        self.assertLess(abs(slope), 0.01)
        self.assertGreater(abs(mismatch[2]+mismatch[0]), 1e-5)

    def test_fixed_rf_parallel_tangent_optimizer(self):
        report = optimize(Settings(532, design_mode="wide-aperture", acoustic_mode="slow",
                                   incident="e", rf_min_mhz=80, rf_max_mhz=80,
                                   rf_budget_w=0.05, grid_deg=5))
        r = report["optimum"]
        self.assertAlmostEqual(r["rf_mhz"], 80, places=6)
        self.assertLess(abs(r["parallel_tangent_residual"]), 1e-10)
        self.assertLess(r["dimensionless_surface_residual"], 1e-12)

    def test_winning_azimuth_response_and_old_report_compatibility(self):
        changed = replace(self.settings, azimuth_deg=37)
        rec = Model(changed).record(0, 15, 0, "o", "e", 1)
        report = dict(inputs=asdict(self.settings), material=asdict(self.material), optimum=rec)
        eta, mismatch = response(report)
        self.assertAlmostEqual(float(eta), rec["throughput_efficiency"], places=10)
        self.assertLess(abs(float(mismatch)), 1e-6)
        old_report = dict(self.report, optimum=dict(self.record))
        old_report["optimum"].pop("azimuth_deg")
        old_eta, _ = response(old_report)
        self.assertAlmostEqual(float(old_eta), self.record["throughput_efficiency"], places=10)


class PackageSearchChecks(unittest.TestCase):
    def cli(self, *args):
        script = Path(__file__).resolve().with_name("teo2_optimizer.py")
        return subprocess.run([sys.executable, str(script), *args],
                              capture_output=True, text=True, timeout=120)

    def test_package_cli_runtime_inputs_and_json(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = str(Path(directory)/"example")
            result = self.cli("--design-package", "--wavelength-nm", "632.8",
                              "--rf-mhz", "60", "--fixed-azimuth",
                              "--design-mode", "wide-aperture", "--acoustic-mode", "slow",
                              "--incident", "e", "--grid-deg", "5", "--output-prefix", prefix)
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(Path(prefix+".json").read_text(),
                                parse_constant=lambda x: self.fail(f"Nonfinite JSON: {x}"))
            self.assertEqual(report["inputs"]["wavelength_nm"], 632.8)
            self.assertAlmostEqual(report["optimum"]["rf_mhz"], 60, places=6)
            manifest = json.loads(Path(prefix+"_design.json").read_text(),
                                  parse_constant=lambda x: self.fail(f"Nonfinite JSON: {x}"))
            self.assertIn("design_status", manifest)
            self.assertEqual(manifest["rf"]["return_loss_target_db"], 30)
            self.assertEqual(manifest["rf"]["reflected_power_fraction_limit"], 0.001)
            self.assertTrue(manifest["rf"]["return_loss_requirement_confirmed"])
            self.assertIn("NOT FOR FABRICATION", result.stdout)
            self.assertTrue(Path(prefix+"_crystal.svg").is_file())
            self.assertTrue(Path(prefix+"_assembly.svg").is_file())

    def test_package_cli_rejects_missing_invalid_conflicting_requirements(self):
        for extra in ([], ["--rf-mhz", "nan"],
                      ["--rf-mhz", "80", "--rf-ceiling-w", "1"],
                      ["--rf-mhz", "80", "--duty-cycle-percent", "0"],
                      ["--rf-mhz", "80", "--duty-cycle-percent", "0.5"],
                      ["--rf-mhz", "80", "--vswr-target", "0.001"],
                      ["--rf-mhz", "80", "--return-loss-target-db", "nan"],
                      ["--rf-mhz", "80", "--measured-vswr", "1000"],
                      ["--rf-mhz", "80", "--rf-min-mhz", "20"],
                      ["--rf-mhz", "80", "--rf-budget-w", "1"]):
            with self.subTest(extra=extra):
                result = self.cli("--design-package", "--wavelength-nm", "532", *extra)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("error:", result.stderr)
        missing_wavelength = self.cli("--design-package", "--rf-mhz", "80")
        self.assertNotEqual(missing_wavelength.returncode, 0)

    def test_bounded_azimuth_fixed_rf_and_feasible_alternatives(self):
        from teo2_design import DesignSettings, crystal_geometry
        engineering = DesignSettings(azimuth_min_deg=40, azimuth_max_deg=50,
                                     azimuth_step_deg=10)
        settings = Settings(532, design_mode="wide-aperture", acoustic_mode="slow",
                            incident="e", rf_min_mhz=80, rf_max_mhz=80, grid_deg=5,
                            length_mm=3, acoustic_width_mm=7, acoustic_height_mm=3)
        report = optimize_design(settings, engineering)
        for rec in [report["optimum"]]+report["alternatives"]:
            self.assertLessEqual(abs(rec["rf_mhz"]-80), 1e-6)
            self.assertTrue(40 <= rec["azimuth_deg"] <= 50)
            candidate = dict(report, optimum=rec)
            self.assertTrue(crystal_geometry(candidate, engineering)["feasible"])
        self.assertEqual(report["inputs"]["azimuth_deg"], report["optimum"]["azimuth_deg"])
        eta, mismatch = response(report)
        self.assertAlmostEqual(float(eta), report["optimum"]["throughput_efficiency"], places=9)
        self.assertLess(abs(float(mismatch)), 1e-5)
        self.assertIn("COPLANAR", report["search"]["domain"])
        self.assertIn(45, report["search"]["azimuth_seeds_deg"])
        self.assertTrue(all(40 <= row["azimuth_deg"] <= 50 for row in report["coarse_cut_scan"]))

    def test_fixed_azimuth_compatibility_and_reject_range(self):
        from teo2_design import DesignSettings
        engineering = DesignSettings(fixed_azimuth=True)
        with self.assertRaisesRegex(ValueError, "fixed nominal"):
            optimize_design(Settings(532), engineering)
        with self.assertRaisesRegex(ValueError, "azimuth bounds"):
            optimize_design(Settings(532, azimuth_deg=100, rf_min_mhz=80, rf_max_mhz=80),
                            engineering)
        settings = Settings(532, design_mode="wide-aperture", acoustic_mode="slow",
                            incident="e", rf_min_mhz=80, rf_max_mhz=80, grid_deg=5,
                            length_mm=3, acoustic_width_mm=7, acoustic_height_mm=3)
        report = optimize_design(settings, engineering)
        self.assertEqual(report["optimum"]["azimuth_deg"], 45)
        self.assertEqual(report["search"]["azimuth_seeds_deg"], [45])

    def test_efficiency_package_search_retains_fixed_rf(self):
        from teo2_design import DesignSettings, crystal_geometry
        engineering = DesignSettings(azimuth_min_deg=40, azimuth_max_deg=50,
                                     azimuth_step_deg=10, rf_ceiling_w=0.2)
        settings = Settings(532, acoustic_mode="slow", incident="e",
                            rf_min_mhz=80, rf_max_mhz=80, grid_deg=10,
                            length_mm=3, acoustic_width_mm=7, acoustic_height_mm=3)
        report = optimize_design(settings, engineering)
        self.assertEqual(report["inputs"]["rf_budget_w"], 0.2)
        for rec in [report["optimum"]]+report["alternatives"]:
            self.assertAlmostEqual(rec["rf_mhz"], 80, places=6)
            self.assertLessEqual(rec["rf_selected_w"], 0.2)
            self.assertTrue(crystal_geometry(dict(report, optimum=rec), engineering)["feasible"])


if __name__ == "__main__":
    unittest.main()
