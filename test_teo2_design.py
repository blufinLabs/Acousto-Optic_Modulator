"""Focused engineering-package checks; these do not qualify a physical device."""
import copy
from dataclasses import asdict, replace
import json
import math
from pathlib import Path
import re
import shutil
import tempfile
import unittest
import xml.etree.ElementTree as ET

import numpy as np

from teo2_design import DesignSettings, air_exit, crystal_geometry, write_design_package
from teo2_optimizer import Material, Model, Settings, response


class DesignChecks(unittest.TestCase):
    def setUp(self):
        self.settings = Settings(632.8, rf_min_mhz=0.01, rf_max_mhz=1000,
                                 rf_budget_w=0.01, optical_transmission=0.95,
                                 length_mm=3, acoustic_width_mm=5, acoustic_height_mm=4)
        self.model = Model(self.settings)
        self.record = self.model.record(0, 15, 0, "o", "e", 1)
        self.assertIsNotNone(self.record)
        self.report = dict(inputs=asdict(self.settings), material=asdict(Material()),
                           optimum=self.record)
        self.engineering = DesignSettings()
        self.artifact_directory = None

    def tearDown(self):
        if self.artifact_directory is not None:
            shutil.rmtree(self.artifact_directory)

    def package(self, engineering=None):
        self.artifact_directory = Path(tempfile.mkdtemp(prefix="teo2-design-", dir="/tmp"))
        return write_design_package(self.report, engineering or self.engineering,
                                    self.artifact_directory/"sample")

    def test_defaults_and_validation(self):
        self.assertEqual(self.engineering.rf_ceiling_w, 0.9)
        self.assertEqual(self.engineering.crystal_length_mm, 10)
        self.assertEqual(self.engineering.vswr_target, 1.5)
        for change in (dict(rf_ceiling_w=0), dict(duty_cycle_percent=101), dict(duty_cycle_percent=0.5),
                       dict(beam_diameter_mm=-1), dict(vswr_target=0.001),
                       dict(measured_vswr=0.1), dict(gaussian_tail_fraction=0),
                       dict(azimuth_step_deg=0), dict(azimuth_max_deg=-1),
                       dict(crystal_width_mm=float("nan")), dict(rf_ceiling_w=True),
                       dict(fixed_azimuth=1), dict(dimension_tolerance_mm=6)):
            with self.subTest(change=change), self.assertRaises(ValueError):
                DesignSettings(**change)
        DesignSettings(measured_vswr=1, azimuth_min_deg=45, azimuth_max_deg=45,
                       fixed_azimuth=True)

    def test_closed_convex_solid_and_unit_normals(self):
        g = crystal_geometry(self.report, self.engineering)
        self.assertTrue(g["feasible"], g["failures"])
        vertices = np.array(g["vertex_coordinates_local_mm"])
        normals = np.array(g["halfspace_normals_local"])
        offsets = np.array(g["halfspace_offsets_mm"])
        self.assertTrue(np.all(normals@vertices.T <= offsets[:, None]+1e-8))
        np.testing.assert_allclose(np.linalg.norm(normals, axis=1), 1, atol=1e-12)
        counts = {}
        for face in g["faces"]:
            n = np.array(face["outward_normal_local"])
            points = vertices[face["vertex_indices"]]
            np.testing.assert_allclose(points@n, face["plane_offset_mm"], atol=1e-8)
            self.assertGreater(face["area_mm2"], 0)
            for a, b in zip(face["vertex_indices"],
                            face["vertex_indices"][1:]+face["vertex_indices"][:1]):
                edge = tuple(sorted((a, b)))
                counts[edge] = counts.get(edge, 0)+1
        self.assertTrue(all(count == 2 for count in counts.values()))
        self.assertEqual(len(vertices)-len(counts)+len(g["faces"]), 2)
        self.assertGreater(g["volume_mm3"], 0)
        self.assertLess(g["volume_mm3"], 1000)
        faces = {f["name"]: f for f in g["faces"]}
        np.testing.assert_allclose(faces["entrance"]["outward_normal_local"], [-1, 0, 0])
        np.testing.assert_allclose(faces["exit"]["outward_normal_local"], [1, 0, 0])

    def test_transforms_and_directed_crystal_frame(self):
        g = crystal_geometry(self.report, asdict(self.engineering))
        a = np.array(g["local_to_crystal_4x4"])
        b = np.array(g["crystal_to_local_4x4"])
        np.testing.assert_allclose(a@b, np.eye(4), atol=1e-12)
        self.assertAlmostEqual(np.linalg.det(a[:3, :3]), 1)
        np.testing.assert_allclose(a[:3, 0], g["incident_wavevector_unit_crystal"])
        np.testing.assert_allclose(np.array(g["vertex_coordinates_local_mm"])@a[:3, :3].T,
                                   g["vertex_coordinates_crystal_mm"])
        for name in ("incident_ray_unit_crystal", "diffracted_ray_unit_crystal",
                     "acoustic_wave_normal_unit_crystal"):
            self.assertAlmostEqual(np.linalg.norm(g[name]), 1)
        self.assertGreater(np.linalg.norm(np.array(g["diffracted_ray_unit_crystal"])-
                                         g["diffracted_wavevector_unit_crystal"]), 0.001)

    def test_alternative_uses_its_own_azimuth_not_winning_inputs(self):
        alternate = Model(replace(self.settings, azimuth_deg=50))
        record = alternate.record(0, 15, 0, "o", "e", 1)
        self.assertIsNotNone(record)
        geometry = crystal_geometry(dict(self.report, optimum=record), self.engineering)
        np.testing.assert_allclose(geometry["incident_wavevector_unit_crystal"],
                                   alternate.direction(15)[0], atol=1e-12)
        old = copy.deepcopy(self.report)
        old["optimum"].pop("azimuth_deg")
        original = crystal_geometry(self.report, self.engineering)
        np.testing.assert_allclose(crystal_geometry(old, self.engineering)["local_to_crystal_4x4"],
                                   original["local_to_crystal_4x4"])

    def test_finite_gaussian_envelope_and_face_footprints(self):
        g = crystal_geometry(self.report, self.engineering)
        beam = g["beam"]
        radius = beam["bounded_radius_mm"]
        self.assertAlmostEqual(math.exp(-2*(radius/(self.engineering.beam_diameter_mm/2))**2),
                               self.engineering.gaussian_tail_fraction)
        self.assertTrue(beam["finite_envelope_only"])
        self.assertTrue(any("diffracted exit" in f["name"] for f in g["optical_footprints"]))
        self.assertTrue(all(f["minimum_clearance_mm"] >= -1e-8 for f in g["optical_footprints"]))
        bad = crystal_geometry(self.report, replace(self.engineering, beam_diameter_mm=4))
        self.assertFalse(bad["feasible"])
        self.assertTrue(any("bounded optical footprint" in f for f in bad["failures"]))

    def test_acoustic_mount_rectangle_on_actual_face_and_group_coverage(self):
        g = crystal_geometry(self.report, self.engineering)
        mount = g["acoustic_mount"]
        corners = np.array(mount["corners_local_mm"])
        normal = np.array(mount["inward_normal_local"])
        center = np.array(mount["center_local_mm"])
        np.testing.assert_allclose((corners-center)@normal, 0, atol=1e-12)
        self.assertAlmostEqual(np.linalg.norm(corners[0]-corners[1]),
                               self.settings.acoustic_height_mm)
        self.assertAlmostEqual(np.linalg.norm(corners[0]-corners[2]),
                               self.settings.acoustic_width_mm)
        self.assertTrue(all(min(c["wavefront_clearances_mm"]) >= 0
                            for c in mount["coverage_checks"]))
        bad = copy.deepcopy(self.report)
        bad["inputs"]["acoustic_height_mm"] = 1
        geometry = crystal_geometry(bad, self.engineering)
        self.assertFalse(geometry["feasible"])
        self.assertTrue(any("bounded acoustic field" in f for f in geometry["failures"]))
        bad["inputs"]["acoustic_width_mm"] = 100
        self.assertTrue(any("mounting" in f for f in
                            crystal_geometry(bad, self.engineering)["failures"]))

    def test_interaction_length_failure(self):
        report = copy.deepcopy(self.report)
        report["inputs"]["length_mm"] = 11
        g = crystal_geometry(report, self.engineering)
        self.assertFalse(g["feasible"])
        self.assertIn("Interaction thickness exceeds optical-face separation", g["failures"])

    def test_fixed_order_finite_continuous_feasibility_margins(self):
        reference = crystal_geometry(self.report, self.engineering)
        names = reference["feasibility_margin_names"]
        values = reference["feasibility_margins"]
        self.assertEqual(len(names), len(values))
        self.assertEqual(len(names), len(set(names)))
        self.assertTrue(all(math.isfinite(x) and x >= 0 for x in values))
        changes = []
        bad = copy.deepcopy(self.report)
        bad["inputs"]["length_mm"] = 11
        changes.append((bad, self.engineering))
        bad = copy.deepcopy(self.report)
        bad["inputs"]["acoustic_height_mm"] = 1
        changes.append((bad, self.engineering))
        bad = copy.deepcopy(self.report)
        bad["inputs"]["acoustic_width_mm"] = 100
        changes.append((bad, self.engineering))
        bad = copy.deepcopy(self.report)
        bad["optimum"]["acoustic_group_velocity_xyz_m_s"] = (
            -np.array(bad["optimum"]["acoustic_normal_xyz"])).tolist()
        changes.append((bad, self.engineering))
        changes.append((self.report, replace(self.engineering, beam_diameter_mm=4)))
        for report, engineering in changes:
            with self.subTest(inputs=report["inputs"], beam=engineering.beam_diameter_mm):
                result = crystal_geometry(report, engineering)
                self.assertFalse(result["feasible"])
                self.assertEqual(result["feasibility_margin_names"], names)
                self.assertTrue(all(math.isfinite(x) for x in result["feasibility_margins"]))
                self.assertTrue(any(x < 0 for x in result["feasibility_margins"]))
        nearby = copy.deepcopy(self.report)
        nearby["optimum"]["incident_wavevector_signed_to_c_deg"] += 1e-5
        result = crystal_geometry(nearby, self.engineering)
        self.assertEqual(result["feasibility_margin_names"], names)
        self.assertLess(np.max(abs(np.array(result["feasibility_margins"])-values)), 1e-4)

    def test_exact_tangential_exit_and_tir(self):
        n = np.array([0., 0., 1.])
        result = air_exit([0.6, 0.2, 2.3], n)
        air = np.array(result["air_direction_xyz"])
        np.testing.assert_allclose(air[:2], [0.6, 0.2])
        self.assertAlmostEqual(air[2], math.sqrt(0.6))
        self.assertAlmostEqual(np.linalg.norm(air), 1)
        self.assertFalse(result["total_internal_reflection"])
        result = air_exit([1.01, 0, 2.3], n)
        self.assertTrue(result["total_internal_reflection"])
        self.assertIsNone(result["air_direction_xyz"])
        self.assertTrue(air_exit([math.sqrt(1+1e-13), 0, 2.3], n)["total_internal_reflection"])
        with self.assertRaises(ValueError):
            air_exit([float("nan"), 0, 0], n)

    def test_anisotropic_output_preserves_crystal_tangential_k(self):
        g = crystal_geometry(self.report, self.engineering)
        n = np.array(g["incident_wavevector_unit_crystal"])
        b = self.record["nd"]*np.array(g["diffracted_wavevector_unit_crystal"])
        result = g["air_outputs"]["diffracted"]
        np.testing.assert_allclose(result["tangential_dimensionless_k"], b-n*(b@n),
                                   atol=1e-12)
        a = np.array(result["air_direction_xyz"])
        np.testing.assert_allclose(a-n*(a@n), b-n*(b@n), atol=1e-12)
        self.assertIsNone(result["power_transmission"])

    def test_six_artifacts_strict_json_and_svg_numeric_projections(self):
        paths = self.package()
        self.assertEqual(len(paths), 6)
        for path in paths.values():
            self.assertTrue(Path(path).exists())
        for name in ("design_json", "crystal_spec_json"):
            text = Path(paths[name]).read_text()
            value = json.loads(text, parse_constant=lambda x: self.fail(f"Nonfinite JSON {x}"))
            self.assertIn("NOT FOR FABRICATION", value["status"])
        for name in ("crystal_drawing", "assembly_schematic"):
            ET.fromstring(Path(paths[name]).read_text())
        text = Path(paths["crystal_drawing"]).read_text()
        self.assertIn("10 mm", text)
        self.assertIn("points=", text)
        self.assertIn("actual clipped solid", text)
        self.assertIn("SCHEMATIC", Path(paths["assembly_schematic"]).read_text())
        self.assertIn("NOT TO SCALE", Path(paths["assembly_schematic"]).read_text())

    def test_fixed_power_rf_tolerance_and_explicit_sensitivity(self):
        paths = self.package()
        data = json.loads(Path(paths["design_json"]).read_text())
        performance = data["modeled_performance"]
        power = self.record["rf_selected_w"]
        self.assertEqual(performance["fixed_rf_power_w"], power)
        rf = performance["rf_carrier_plus_minus_0p1_percent"]
        self.assertEqual([x["fractional_offset"] for x in rf], [-0.001, 0, 0.001])
        for point in rf:
            eta, mismatch = response(self.report, rf_mhz=point["rf_mhz"], rf_power_w=power)
            self.assertAlmostEqual(point["throughput_efficiency_fraction"], float(eta))
            self.assertAlmostEqual(point["longitudinal_mismatch_m_inv"], float(mismatch))
        self.assertEqual(len(performance["acoustic_cut"]), 2)
        self.assertEqual(len(performance["internal_alignment"]), 3)
        self.assertEqual(len(performance["crystal_dimensions"]), 6)
        self.assertGreater(abs(performance["acoustic_cut"][0]["longitudinal_mismatch_m_inv"]), 0.1)

    def test_no_false_release_or_invented_match_or_performance(self):
        paths = self.package()
        data = json.loads(Path(paths["design_json"]).read_text())
        self.assertFalse(data["fabrication_release"])
        self.assertIsNone(data["rf"]["net_accepted_rf_w"])
        self.assertIsNone(data["rf"]["measured_vswr"])
        self.assertTrue(data["rf"]["requirement_history"]["superseded"])
        self.assertFalse(data["rf"]["requirement_history"]["reversal_valid"])
        self.assertTrue(data["rf"]["return_loss_requirement_confirmed"])
        self.assertEqual(data["rf"]["requirement_history"]["subsequent"], "1:1000")
        self.assertTrue(all(x is None for x in data["unmodeled_performance"].values()))
        self.assertTrue(all(x["evidence"] is None for x in data["acceptance"]))
        self.assertGreaterEqual(len(data["release_blockers"]), 10)
        self.assertTrue(any("LiNbO3" in x["statement"] for x in data["release_blockers"]))
        self.assertTrue(any("Thermal" in x["statement"] for x in data["release_blockers"]))
        self.assertIn("20-10", Path(paths["design_report"]).read_text())

    def test_tolerance_proposals_have_rationale_and_no_fabricated_source_evidence(self):
        paths = self.package()
        data = json.loads(Path(paths["design_json"]).read_text())
        tolerances = data["proposed_tolerances"]
        proposals = [tolerances[name] for name in ("dimension", "acoustic_cut",
                                                   "internal_alignment", "substrate_orientation",
                                                   "optical_face_orientation")]
        proposals += list(tolerances["surface_quality"]["per_value_proposals"].values())
        for item in proposals:
            self.assertIn("UNVERIFIED", item["source_status"])
            self.assertTrue(item["selection_rationale"])
            self.assertTrue(item["attempted_primary_sources"])
        sources = tolerances["primary_source_verification"]["sources_to_verify"]
        self.assertEqual(len(sources), 3)
        self.assertTrue(all(x["url"].startswith("https://") and
                            not x["numerical_claim_verified"] for x in sources))
        self.assertTrue(any("primary-source" in blocker["statement"] for blocker in data["release_blockers"]))
        flatness = tolerances["surface_quality"]["per_value_proposals"]["surface_flatness"]
        self.assertAlmostEqual(flatness["value"]*flatness["reference_wavelength_nm"],
                               flatness["equivalent_surface_height_pv_nm"])
        self.assertEqual(flatness["reference_aperture"]["datum_faces"], ["A", "B"])
        self.assertEqual(len(flatness["reference_aperture"]["footprints"]), 4)
        roughness = tolerances["surface_quality"]["per_value_proposals"]["roughness"]
        self.assertIn("Sq", roughness["metric"])
        for field in ("value", "spatial_bandwidth_mm_inv", "sampling_pitch_um", "scan_size_mm"):
            self.assertIsNone(roughness[field])

    def test_drawing_contains_numeric_datums_axes_apertures_and_cut_metadata(self):
        paths = self.package()
        data = json.loads(Path(paths["design_json"]).read_text())
        document = ET.fromstring(Path(paths["crystal_drawing"]).read_text())
        text = " ".join(document.itertext())
        for value in ("CONCEPT-01", "[100]", "[010]", "[001]", "DATUMS", "FACE PLANES",
                      "TRANSFORM", "CUT ANGLES", "INTERACTION SLAB", "TRANSDUCER ACTIVE AREA",
                      "TOTAL CERAMIC", "BOUNDED APERTURES", "VERTICES", "Roughness"):
            self.assertIn(value, text)
        self.assertGreater(int(document.attrib["height"]), 1000)
        rotation = np.array(data["geometry"]["local_to_crystal_4x4"])[:3, :3]
        axis_paths = [item for item in document.iter() if "data-crystal-axis" in item.attrib]
        self.assertEqual(len(axis_paths), 9)
        for item in axis_paths:
            crystal_axis = {"[100]": 0, "[010]": 1, "[001]": 2}[item.attrib["data-crystal-axis"]]
            a, b = ("xyz".index(x) for x in item.attrib["data-view"])
            numbers = [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", item.attrib["d"])]
            projected = rotation.T@np.eye(3)[crystal_axis]
            np.testing.assert_allclose([numbers[2]-numbers[0], numbers[3]-numbers[1]],
                                       [35*projected[a], -35*projected[b]], atol=1e-5)

    def test_state_blocker_acceptance_and_continuous_duty_are_not_toggles(self):
        paths = self.package(replace(self.engineering, duty_cycle_percent=10))
        data = json.loads(Path(paths["design_json"]).read_text())
        self.assertEqual(data["design_status"], "preliminary_concept")
        self.assertEqual(data["allowed_design_states"],
                         ["preliminary_concept", "review_ready", "fabrication_released"])
        self.assertTrue(data["design_status_reasons"])
        self.assertFalse(data["fabrication_release"])
        self.assertIsNone(data["conditional_assumptions"]["interaction_slab_physically_realized"])
        self.assertTrue(any(b["id"] == "physical-interaction-length" for b in data["release_blockers"]))
        for blocker in data["release_blockers"]:
            self.assertTrue(blocker["verification"])
            self.assertTrue(blocker["acceptance_criterion"])
            self.assertEqual(blocker["status"], "open")
            self.assertIsNone(blocker["evidence"])
        thermal = next(b for b in data["release_blockers"] if b["id"] == "continuous-duty-thermal")
        self.assertIn("100%", thermal["verification"])
        self.assertIn("regardless", thermal["verification"])
        self.assertIn("continuous 100%", next(t for t in data["acceptance"]
                                              if t["test"] == "Thermal / endurance")["criterion"])

    def test_interaction_volume_bound_and_active_area_are_not_ceramic_outline(self):
        g = crystal_geometry(self.report, self.engineering)
        slab = g["interaction_slab"]
        self.assertEqual(slab["thickness_mm"], 3)
        self.assertEqual(slab["crystal_optical_face_separation_mm"], 10)
        lower, upper = np.array(slab["bounded_optical_envelope_box_local_mm"])
        self.assertAlmostEqual(slab["bounded_optical_volume_upper_bound_mm3"],
                               np.prod(upper-lower))
        radius = g["beam"]["bounded_radius_mm"]
        self.assertGreaterEqual(slab["bounded_optical_volume_upper_bound_mm3"],
                                math.pi*radius**2*slab["thickness_mm"])
        self.assertEqual(g["acoustic_mount"]["active_acoustic_area_mm2"], 20)
        self.assertIsNone(g["acoustic_mount"]["total_ceramic_dimensions_mm"])

    def test_measured_match_distinct_from_target(self):
        paths = self.package(replace(self.engineering, measured_vswr=2, vswr_target=1.1))
        data = json.loads(Path(paths["design_json"]).read_text())
        self.assertAlmostEqual(data["rf"]["net_accepted_rf_w"], self.record["rf_selected_w"]*8/9)
        self.assertEqual(data["rf"]["target_vswr"], 1.1)
        self.assertEqual(data["rf"]["measured_vswr"], 2)
        self.assertAlmostEqual(data["rf"]["measured_reflection_magnitude"], 1/3)
        self.assertFalse(data["rf"]["return_loss_meets_target_from_supplied_vswr"])

    def test_rf_match_requirement_is_inclusive_band_limited_and_not_optical_loss(self):
        paths = self.package()
        data = json.loads(Path(paths["design_json"]).read_text())
        rf = data["rf"]
        self.assertEqual(rf["return_loss_target_db"], 30)
        self.assertEqual(rf["reflected_power_fraction_limit"], 0.001)
        self.assertAlmostEqual(rf["equivalent_vswr_limit"], 1.0653108640674351)
        self.assertEqual(rf["stated_vswr_limit"], 1.065)
        self.assertEqual(rf["return_loss_comparison"], "greater than or equal to")
        self.assertEqual(rf["reflected_power_fraction_comparison"], "less than or equal to")
        self.assertEqual(rf["s11_db_upper_limit"], -30)
        self.assertEqual(rf["rf_match_frequency_band"],
                         "[f_min, f_max]; numerical endpoints to be agreed")
        self.assertEqual(rf["port2_termination_ohm"], 50)
        self.assertIsNone(rf["measured_return_loss_db"])
        self.assertIsNone(rf["return_loss_meets_target_from_supplied_vswr"])
        self.assertIsNone(rf["net_accepted_rf_w"])
        self.assertFalse(rf["measured_vswr_is_band_compliance_evidence"])
        self.assertTrue(rf["return_loss_requirement_confirmed"])
        for value in (0, -1, 301, float("nan")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                DesignSettings(return_loss_target_db=value)
        shutil.rmtree(self.artifact_directory)
        self.artifact_directory = None
        paths = self.package(replace(self.engineering, measured_vswr=1))
        data = json.loads(Path(paths["design_json"]).read_text())
        self.assertTrue(data["rf"]["supplied_vswr_is_zero_reflection_limit"])
        self.assertIsNone(data["rf"]["measured_return_loss_db"])  # Infinite ideal limit, not JSON infinity.
        self.assertTrue(data["rf"]["return_loss_meets_target_from_supplied_vswr"])
        self.assertTrue(data["rf"]["measured_vswr_meets_stated_limit"])
        self.assertTrue(data["rf"]["measured_point_meets_rf_limits"])
        self.assertFalse(data["fabrication_release"])
        shutil.rmtree(self.artifact_directory)
        self.artifact_directory = None
        at_limit = (1+math.sqrt(0.001))/(1-math.sqrt(0.001))
        paths = self.package(replace(self.engineering, measured_vswr=at_limit))
        data = json.loads(Path(paths["design_json"]).read_text())
        self.assertTrue(data["rf"]["return_loss_meets_target_from_supplied_vswr"])
        self.assertFalse(data["rf"]["measured_vswr_meets_stated_limit"])
        self.assertFalse(data["rf"]["measured_point_meets_rf_limits"])
        shutil.rmtree(self.artifact_directory)
        self.artifact_directory = None
        paths = self.package(replace(self.engineering, measured_vswr=1.065))
        data = json.loads(Path(paths["design_json"]).read_text())
        self.assertTrue(data["rf"]["measured_point_meets_rf_limits"])
        self.assertFalse(data["rf"]["measured_vswr_is_band_compliance_evidence"])

    def test_ceiling_violation_stays_blocked(self):
        paths = self.package(replace(self.engineering, rf_ceiling_w=0.001))
        data = json.loads(Path(paths["design_json"]).read_text())
        self.assertFalse(data["rf"]["ceiling_met"])
        self.assertIn("Selected RF power exceeds engineering RF ceiling",
                      [b["statement"] for b in data["release_blockers"]])
        self.assertFalse(data["fabrication_release"])


if __name__ == "__main__":
    unittest.main()
