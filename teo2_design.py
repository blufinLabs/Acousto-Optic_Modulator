"""Concept engineering package for the conditional bulk TeO2 model.

Geometry is in millimetres. Crystal Cartesian coordinates are [100], [010],
[001]; drawing coordinates have x along the incident wavevector. No artifact
produced here is a manufacturing release or a substitute for measured evidence.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields, replace
import copy
import itertools
import json
import math
from pathlib import Path
import textwrap
from xml.sax.saxutils import escape

import numpy as np

DESIGN_REVISION = "CONCEPT-01"
DESIGN_STATES = ["preliminary_concept", "review_ready", "fabrication_released"]


@dataclass(frozen=True)
class DesignSettings:
    rf_ceiling_w: float = 0.9
    duty_cycle_percent: float = 100.0
    crystal_length_mm: float = 10.0
    crystal_width_mm: float = 10.0
    crystal_height_mm: float = 10.0
    beam_diameter_mm: float = 1.0
    vswr_target: float = 1.5
    azimuth_min_deg: float = 0.0
    azimuth_max_deg: float = 90.0
    azimuth_step_deg: float = 15.0
    fixed_azimuth: bool = False
    gaussian_tail_fraction: float = 1e-6
    measured_vswr: float = 0.0
    dimension_tolerance_mm: float = 0.05
    cut_tolerance_deg: float = 0.05
    alignment_tolerance_deg: float = 0.01
    return_loss_target_db: float = 30.0

    def __post_init__(self):
        self.validate()

    def validate(self):
        for f in fields(self):
            value = getattr(self, f.name)
            if f.name == "fixed_azimuth":
                if not isinstance(value, bool):
                    raise ValueError("fixed_azimuth must be boolean")
            elif isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"{f.name} must be a finite number")
        for name in ("rf_ceiling_w", "crystal_length_mm", "crystal_width_mm",
                     "crystal_height_mm", "beam_diameter_mm", "azimuth_step_deg",
                     "dimension_tolerance_mm", "cut_tolerance_deg", "alignment_tolerance_deg"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")
        if self.rf_ceiling_w >= 1:
            raise ValueError("rf_ceiling_w must be strictly below 1 W for this subwatt design")
        if not 1 <= self.duty_cycle_percent <= 100:
            raise ValueError("duty_cycle_percent must lie in [1,100]")
        if self.vswr_target < 1:
            raise ValueError("VSWR is Vmax/Vmin and must be >=1; 1:1000 is not valid VSWR")
        if not 0 < self.return_loss_target_db <= 300:
            raise ValueError("return_loss_target_db must lie in (0,300] dB for finite floating-point RF ratios")
        if self.measured_vswr != 0 and self.measured_vswr < 1:
            raise ValueError("measured_vswr must be 0 (not supplied) or >=1")
        if not 0 <= self.azimuth_min_deg <= self.azimuth_max_deg <= 360:
            raise ValueError("Azimuth bounds must increase within [0,360]")
        if not self.fixed_azimuth and self.azimuth_min_deg == self.azimuth_max_deg:
            raise ValueError("An azimuth search requires a nonzero range")
        if not 0 < self.gaussian_tail_fraction < 1:
            raise ValueError("gaussian_tail_fraction must lie in (0,1)")
        if self.dimension_tolerance_mm >= min(self.crystal_length_mm, self.crystal_width_mm,
                                              self.crystal_height_mm)/2:
            raise ValueError("Dimension tolerance exceeds crystal half-size")
        if max(self.cut_tolerance_deg, self.alignment_tolerance_deg) >= 5:
            raise ValueError("Perturbation angles must be below 5 degrees")


def _engineering(value):
    if isinstance(value, DesignSettings):
        value.validate()
        return value
    return DesignSettings(**value)


def _unit(value):
    v = np.asarray(value, dtype=float)
    if v.shape != (3,) or not np.all(np.isfinite(v)) or np.linalg.norm(v) < 1e-12:
        raise ValueError("Expected a finite nonzero three-vector")
    return v/np.linalg.norm(v)


def air_exit(dimensionless_k, outward_normal):
    """Conserve tangential k/k0, not a scalar extraordinary refractive index."""
    b, n = np.asarray(dimensionless_k, dtype=float), _unit(outward_normal)
    if b.shape != (3,) or not np.all(np.isfinite(b)):
        raise ValueError("dimensionless_k must be a finite three-vector")
    bt = b-n*np.dot(b, n)
    t2 = float(bt@bt)
    tir = t2 > 1.0
    k_air = None if tir else bt+n*math.sqrt(max(0.0, 1.0-t2))
    return dict(tangential_dimensionless_k=bt.tolist(), tangential_norm_squared=t2,
                total_internal_reflection=tir,
                air_direction_xyz=None if k_air is None else k_air.tolist(),
                grazing=not tir and 1-t2-2e-12 < 0,
                boundary_equation="b_air = b_t + n*sqrt(1 - |b_t|^2); n_air=1",
                power_transmission=None,
                note="Direction only: Fresnel amplitudes, coating, optical activity and RF optical-frequency shift omitted.")


def _solid(normals, offsets, labels):
    vertices = []
    for ix in itertools.combinations(range(len(normals)), 3):
        a = normals[list(ix)]
        if abs(np.linalg.det(a)) < 1e-10:
            continue
        p = np.linalg.solve(a, offsets[list(ix)])
        if np.all(normals@p <= offsets+1e-8) and not any(
                np.linalg.norm(p-old) < 1e-7 for old in vertices):
            vertices.append(p)
    if len(vertices) < 4:
        raise ValueError("Clipped crystal is empty or degenerate")
    vertices = np.array(vertices)
    faces, edges, volume = [], {}, 0.0
    for i, n in enumerate(normals):
        if any(np.linalg.norm(n-normals[j]) < 1e-10 and
               abs(offsets[i]-offsets[j]) < 1e-10 for j in range(i+1, len(normals))):
            continue
        indices = np.flatnonzero(abs(vertices@n-offsets[i]) < 1e-7)
        if len(indices) < 3:
            continue
        center = np.mean(vertices[indices], axis=0)
        axis = _unit(vertices[indices[0]]-center)
        other = np.cross(n, axis)
        order = np.argsort(np.arctan2((vertices[indices]-center)@other,
                                     (vertices[indices]-center)@axis))
        indices = indices[order].tolist()
        points = vertices[indices]
        area = 0.0
        for j in range(len(points)):
            area += float(np.dot(np.cross(points[j]-center,
                                          points[(j+1) % len(points)]-center), n))/2
            edge = tuple(sorted((indices[j], indices[(j+1) % len(points)])))
            edges[edge] = edges.get(edge, 0)+1
        for j in range(1, len(points)-1):
            volume += float(np.dot(points[0], np.cross(points[j], points[j+1])))/6
        faces.append(dict(name=labels[i], vertex_indices=indices, outward_normal_local=n.tolist(),
                          plane_offset_mm=float(offsets[i]), area_mm2=area,
                          centroid_local_mm=center.tolist()))
    if any(count != 2 for count in edges.values()) or volume <= 1e-9:
        raise ValueError("Crystal faces do not form a closed positive-volume solid")
    return vertices, faces, volume, sorted(edges)


def crystal_geometry(report, engineering):
    """Construct and check a clipped box and bounded optical/acoustic footprints.

    Gaussian diameter means intensity 1/e^2 diameter. Its infinite tails cannot
    be contained: feasibility concerns the explicitly reported power envelope,
    not literal containment of a Gaussian or a calculated overlap integral.
    """
    from teo2_optimizer import Model, Settings
    e = _engineering(engineering)
    s = Settings(**report["inputs"])
    r = report["optimum"]
    s = replace(s, azimuth_deg=r.get("azimuth_deg", s.azimuth_deg))
    model = Model(s)
    theta = r["incident_wavevector_signed_to_c_deg"]
    si = _unit(model.direction(theta)[0])
    sd = _unit(model.direction(r["diffracted_wavevector_signed_to_c_deg"])[0])
    ni = float(r["ni"])
    nd = float(r["nd"])
    if min(ni, nd) <= 0 or not all(math.isfinite(v) for v in (ni, nd)):
        raise ValueError("Optical indices must be finite and positive")
    bi, bd = ni*si, nd*sd
    ri = _unit(model.surfaces[r["incident_polarization"]]@bi)
    rd = _unit(model.surfaces[r["diffracted_polarization"]]@bd)
    z = _unit(model.o)
    y = _unit(np.cross(z, si))
    rotation = np.column_stack((si, y, z))
    q_crystal = _unit(r["acoustic_normal_xyz"])
    q = rotation.T@q_crystal
    if abs(q[2]) > 1e-7:
        raise ValueError("Acoustic normal must lie in the specified interaction plane")
    q[2] = 0.0
    q = _unit(q)
    tangent = _unit(np.cross([0., 0., 1.], q))
    vg = _unit(rotation.T@np.asarray(r["acoustic_group_velocity_xyz_m_s"], dtype=float))
    advance = float(q@vg)
    # Place the active-area center on the group ray through the interaction
    # center, rather than on the wave normal: these differ at acoustic walkoff.
    mounting_direction = vg/advance if advance > 1e-10 else q
    half = np.array([e.crystal_length_mm, e.crystal_width_mm, e.crystal_height_mm])/2
    width, height = s.acoustic_width_mm, s.acoustic_height_mm
    failures = []
    constraint_names, constraint_values = [], []

    def constraint(name, value, tolerance=0.0):
        value = float(value)
        if not math.isfinite(value):
            value = -1.0
        elif -tolerance <= value < 0:
            value = 0.0
        constraint_names.append(name)
        constraint_values.append(value)

    constraint("interaction_thickness_fits_mm", e.crystal_length_mm-s.length_mm)
    if s.length_mm > e.crystal_length_mm:
        failures.append("Interaction thickness exceeds optical-face separation")
    margins = half-abs(tangent)*width/2-np.array([0., 0., height/2])
    bounds = [margins[j]/abs(mounting_direction[j])
              for j in range(3) if abs(mounting_direction[j]) > 1e-10]
    distance = min(bounds)
    for axis, margin in zip("xyz", margins):
        constraint(f"active_wavefront_box_half_margin_{axis}_mm", margin)
    constraint("mount_plane_inward_distance_mm", distance-1e-12)
    if distance <= 0 or np.min(margins) < 0:
        failures.append("Acoustic wavefront rectangle cannot fit an oblique mounting facet")
    distance = max(1e-4, distance)
    normals = np.array([[-1., 0, 0], [1., 0, 0], [0, -1., 0],
                        [0, 1., 0], [0, 0, -1.], [0, 0, 1.], -q])
    offsets = np.r_[half[0], half[0], half[1], half[1], half[2], half[2], distance]
    labels = ["entrance", "exit", "side_minus_y", "side_plus_y",
              "side_minus_z", "side_plus_z", "acoustic_mount"]
    vertices, faces, volume, edges = _solid(normals, offsets, labels)
    for face in faces:
        face["outward_normal_crystal"] = (rotation@face["outward_normal_local"]).tolist()
        face["centroid_crystal_mm"] = (rotation@face["centroid_local_mm"]).tolist()
    mount_center = -mounting_direction*distance
    corners = np.array([mount_center+a*width/2*tangent+b*height/2*np.array([0., 0, 1.])
                        for a, b in itertools.product((-1, 1), repeat=2)])
    mount_clearance = float(np.min(offsets[:, None]-normals@corners.T))
    mount_face = next((f for f in faces if f["name"] == "acoustic_mount"), None)
    if mount_clearance < -1e-7 or mount_face is None:
        failures.append("Specified acoustic wavefront extends outside the mounting polygon")
    if abs(q[0]) > 1-1e-8:
        failures.append("Acoustic mount coincides with an optical face")
    constraint("active_rectangle_mount_corner_clearance_mm", mount_clearance, 1e-7)
    constraint("mount_facet_area_mm2", -1.0 if mount_face is None else mount_face["area_mm2"])
    constraint("acoustic_mount_distinct_from_optical_faces", 1-abs(q[0])-1e-8)
    radius = e.beam_diameter_mm/2*math.sqrt(-math.log(e.gaussian_tail_fraction)/2)
    local_ri, local_rd = rotation.T@ri, rotation.T@rd
    if min(local_ri[0], local_rd[0]) <= 1e-10:
        failures.append("Optical energy ray is not forward through the parallel optical faces")
    constraint("incident_ray_forward_projection", local_ri[0]-1e-10)
    constraint("diffracted_ray_forward_projection", local_rd[0]-1e-10)
    constraint("acoustic_group_inward_projection", advance-1e-10)
    footprints, paths, interaction_samples = [], [], []
    min_margin = math.inf

    def section(ray, x, origin):
        if abs(ray[0]) < 1e-10:
            raise ValueError("Optical ray is tangent to an optical face")
        center = origin+ray*((x-origin[0])/ray[0])
        projection = np.eye(3)-np.outer(ray, [1., 0., 0.])/ray[0]
        projector = np.eye(3)-np.outer(ray, ray)
        spread = radius*np.linalg.norm(projector@projection.T@normals.T, axis=0)
        margin = offsets-normals@center-spread
        return center, projection, projector, margin

    def check(name, ray, x, origin, interaction=False):
        nonlocal min_margin
        center, projection, projector, margin = section(ray, x, origin)
        minimum = float(np.min(margin))
        # A footprint on its own optical face has identically zero normal
        # clearance; report usable edge margin instead of that tautology.
        applicable = [i for i in range(len(normals)) if not (
            i < 2 and abs(center[0]-(-half[0] if i == 0 else half[0])) < 1e-8)]
        edge_margin = float(np.min(margin[applicable]))
        min_margin = min(min_margin, edge_margin)
        if minimum < -1e-7:
            failures.append(f"{name} bounded optical footprint crosses a crystal edge ({minimum:.6g} mm)")
        constraint(name+"_optical_edge_clearance_mm", np.min(margin[2:]), 1e-7)
        item = dict(name=name, center_local_mm=center.tolist(), minimum_clearance_mm=minimum,
                    clear_aperture_edge_margin_mm=edge_margin,
                    envelope_projection_local=projection.tolist(),
                    transverse_projector_local=projector.tolist())
        transverse_axis = _unit(np.cross(ray, [0., 0., 1.]))
        item["envelope_axis1_local_mm"] = (radius*projection@transverse_axis).tolist()
        item["envelope_axis2_local_mm"] = (radius*projection@np.cross(ray, transverse_axis)).tolist()
        item["envelope_half_extents_local_mm"] = (
            radius*np.linalg.norm(projector@projection.T, axis=0)).tolist()
        footprints.append(item)
        if interaction:
            interaction_samples.append((name, center, projection, projector))
        return center

    origin = np.zeros(3)
    entrance = check("incident entrance", local_ri, -half[0], origin)
    zero_exit = check("zero-order exit", local_ri, half[0], origin)
    paths.append(dict(name="incident / zero order", start_local_mm=entrance.tolist(),
                      end_local_mm=zero_exit.tolist(), ray_local=local_ri.tolist()))
    for sign in (-1, 1):
        x = sign*s.length_mm/2
        source = check(f"incident interaction {sign:+d}", local_ri, x, origin, True)
        end = check(f"diffracted exit from interaction {sign:+d}", local_rd, half[0], source)
        paths.append(dict(name=f"diffracted from slab {sign:+d}", start_local_mm=source.tolist(),
                          end_local_mm=end.tolist(), ray_local=local_rd.tolist()))
        check(f"diffracted interaction {sign:+d}", local_rd, s.length_mm/2, source, True)
    acoustic_checks = []
    if advance <= 1e-10:
        failures.append("Acoustic group velocity does not propagate inward from mount")
        for name, *_ in interaction_samples:
            for suffix in ("acoustic_width_clearance_mm", "acoustic_height_clearance_mm",
                           "acoustic_forward_margin_mm"):
                constraint(name+"_"+suffix, -1.0)
    else:
        launch = np.eye(3)-np.outer(vg, q)/advance
        for name, center, projection, projector in interaction_samples:
            foot = center-vg*(q@(center-mount_center)/advance)
            axes = np.array([tangent, [0., 0., 1.]])
            extent = radius*np.linalg.norm(projector@projection.T@launch.T@axes.T, axis=0)
            coordinate = axes@(foot-mount_center)
            clearance = np.array([width, height])/2-abs(coordinate)-extent
            forward_margin = float(q@(center-mount_center)-radius*np.linalg.norm(
                projector@projection.T@q))
            minimum = min(float(np.min(clearance)), forward_margin)
            acoustic_checks.append(dict(name=name, launch_center_local_mm=foot.tolist(),
                                        wavefront_clearances_mm=clearance.tolist(),
                                        forward_margin_mm=forward_margin))
            if minimum < -1e-7:
                failures.append(f"{name} leaves bounded acoustic field including group walkoff")
            constraint(name+"_acoustic_width_clearance_mm", clearance[0], 1e-7)
            constraint(name+"_acoustic_height_clearance_mm", clearance[1], 1e-7)
            constraint(name+"_acoustic_forward_margin_mm", forward_margin, 1e-7)
    exits = {"zero_order": air_exit(bi, si), "diffracted": air_exit(bd, si)}
    for name, result in exits.items():
        if result["total_internal_reflection"] or result["grazing"]:
            failures.append(f"{name} output is totally internally reflected or grazing in air")
        constraint(name+"_air_normal_wavevector_squared",
                   1-result["tangential_norm_squared"]-2e-12)
    transform = np.eye(4)
    transform[:3, :3] = rotation
    inverse = transform.T.copy()
    lower, upper = np.full(3, math.inf), np.full(3, -math.inf)
    for _, center, projection, projector in interaction_samples:
        extent = radius*np.linalg.norm(projector@projection.T, axis=0)
        lower = np.minimum(lower, center-extent)
        upper = np.maximum(upper, center+extent)
    lower[0], upper[0] = -s.length_mm/2, s.length_mm/2
    slab = dict(boundary_normal_local=[1., 0., 0.], thickness_mm=s.length_mm,
                crystal_optical_face_separation_mm=e.crystal_length_mm,
                x_interval_mm=[-s.length_mm/2, s.length_mm/2],
                bounded_optical_envelope_box_local_mm=[lower.tolist(), upper.tolist()],
                bounded_optical_volume_upper_bound_mm3=float(np.prod(upper-lower)),
                volume_definition="Conservative ray-envelope bounding prism over the interaction slab, not actual Gaussian intensity volume or an optical/acoustic overlap integral. Its corners may extend outside the clipped crystal; feasibility tests the ray envelopes, not this prism.",
                volume_validity="Valid for the explicitly bounded constant-radius ray model only; not a diffraction/divergence calculation.",
                physical_boundary_status="Model slab boundaries only, not a verified acoustic-field cutoff. The active wavefront can cause coupling outside this slab; physical localization/effective interaction length is a release blocker.")
    return dict(feasible=not failures, failures=list(dict.fromkeys(failures)),
                feasibility_margins=constraint_values,
                feasibility_margin_names=constraint_names,
                feasibility_margin_convention="Fixed-order finite residuals, >=0 feasible; mm unless names specify dimensionless projection, wavevector squared or mm2. Roundoff within the containment tolerance is clamped to zero. Min/absolute-value bounds are continuous but piecewise differentiable.",
                units="mm", construction="intersection of six nominal box halfspaces and oblique acoustic-mount halfspace",
                nominal_dimensions_mm=(half*2).tolist(), vertex_coordinates_local_mm=vertices.tolist(),
                vertex_coordinates_crystal_mm=(vertices@rotation.T).tolist(),
                faces=faces, edges=[list(x) for x in edges], volume_mm3=volume,
                halfspace_normals_local=normals.tolist(), halfspace_offsets_mm=offsets.tolist(),
                local_to_crystal_4x4=transform.tolist(), crystal_to_local_4x4=inverse.tolist(),
                transform_validation=dict(orthogonality_residual=float(np.max(abs(rotation.T@rotation-np.eye(3)))),
                                          determinant=float(np.linalg.det(rotation))),
                incident_wavevector_unit_crystal=si.tolist(), diffracted_wavevector_unit_crystal=sd.tolist(),
                incident_ray_unit_crystal=ri.tolist(), diffracted_ray_unit_crystal=rd.tolist(),
                acoustic_wave_normal_unit_crystal=q_crystal.tolist(), acoustic_group_unit_local=vg.tolist(),
                cut_angles=dict(azimuth_deg=s.azimuth_deg,
                                acoustic_normal_to_plus_c_deg=math.degrees(math.acos(float(np.clip(q_crystal[2], -1, 1)))),
                                acoustic_face_plane_to_c_acute_deg=math.degrees(math.asin(abs(float(q_crystal[2])))),
                                incident_wavevector_signed_to_c_deg=theta,
                                incident_wavevector_to_plus_c_deg=math.degrees(math.acos(float(np.clip(si[2], -1, 1)))),
                                diffracted_wavevector_signed_to_c_deg=r["diffracted_wavevector_signed_to_c_deg"]),
                optical_polarizations=dict(incident=r["incident_polarization"],
                                           diffracted=r["diffracted_polarization"]),
                operating_point=dict(wavelength_nm=s.wavelength_nm, rf_mhz=r["rf_mhz"],
                                     diffraction_frequency_shift_sign=r["frequency_shift_sign"],
                                     acoustic_mode_index=r["acoustic_mode_index"],
                                     return_loss_target_db=e.return_loss_target_db),
                optical_faces_parallel=True, entrance_normal_parallel_incident_k=True,
                beam=dict(diameter_definition="intensity 1/e^2 diameter", diameter_mm=e.beam_diameter_mm,
                          bounded_radius_mm=radius, excluded_power_fraction=e.gaussian_tail_fraction,
                          contained_power_fraction=1-e.gaussian_tail_fraction,
                          finite_envelope_only=True,
                          caveat="Gaussian tails are infinite. Diameter is internal in the plane perpendicular to each energy ray. Constant-radius collimated ray envelope; diffraction/divergence, astigmatism and integrated overlap not modeled."),
                optical_footprints=footprints, optical_paths=paths,
                interaction_slab=slab,
                minimum_optical_clearance_mm=min_margin,
                acoustic_mount=dict(center_local_mm=mount_center.tolist(), inward_normal_local=q.tolist(),
                                    width_axis_local=tangent.tolist(), height_axis_local=[0., 0., 1.],
                                    wavefront_dimensions_mm=[width, height], corners_local_mm=corners.tolist(),
                                    active_acoustic_area_mm2=width*height,
                                    total_ceramic_dimensions_mm=None, ceramic_thickness_mm=None,
                                    electrode_outline_mm=None, ceramic_border_mm=None,
                                    note="Fitted rectangle is the modeled acoustic active area only. Total ceramic footprint, electrode border, bond and interconnect clearance are unresolved; active-area fit is not ceramic fit.",
                                    minimum_corner_clearance_mm=mount_clearance,
                                    acoustic_field_model="uniform rectangular wavefront advected along acoustic group velocity",
                                    coverage_checks=acoustic_checks),
                air_outputs=exits)


def _finite(value):
    """Represent unavailable model evaluations as null, never JSON NaN."""
    if isinstance(value, np.ndarray):
        return _finite(value.tolist())
    if isinstance(value, dict):
        return {str(k): _finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def _performance(report, e):
    from teo2_optimizer import Material, acoustic_modes, response
    r, s = report["optimum"], report["inputs"]
    power, frequency = r["rf_selected_w"], r["rf_mhz"]

    def evaluate(data, **kwargs):
        eta, mismatch = response(data, rf_power_w=power, **kwargs)
        eta, mismatch = float(eta), float(mismatch)
        return dict(throughput_efficiency_fraction=eta, longitudinal_mismatch_m_inv=mismatch,
                    evaluation_valid=math.isfinite(eta) and math.isfinite(mismatch))

    rf = [dict(fractional_offset=x, rf_mhz=frequency*(1+x),
               **evaluate(report, rf_mhz=frequency*(1+x))) for x in (-0.001, 0, 0.001)]
    alignment = [dict(internal_wavevector_offset_deg=x,
                      **evaluate(report, alignment_offset_deg=x))
                 for x in (-e.alignment_tolerance_deg, 0, e.alignment_tolerance_deg)]
    cut = []
    phi = math.radians(r.get("azimuth_deg", s["azimuth_deg"]))
    basal = np.array([math.cos(phi), math.sin(phi), 0.])
    previous_u = _unit(r["acoustic_displacement_xyz"])
    material = Material(**report["material"])
    for delta in (-e.cut_tolerance_deg, e.cut_tolerance_deg):
        alpha = math.radians(r["acoustic_tilt_from_basal_deg"]+delta)
        q = math.cos(alpha)*basal+math.sin(alpha)*np.array([0., 0., 1.])
        speeds, eigenvectors = acoustic_modes(q, material)
        mode = int(np.argmax(abs(eigenvectors.T@previous_u)))
        u = eigenvectors[:, mode]
        if u@previous_u < 0:
            u = -u
        changed = copy.deepcopy(report)
        changed["optimum"].update(acoustic_normal_xyz=q.tolist(), acoustic_displacement_xyz=u.tolist(),
                                  acoustic_phase_speed_m_s=float(speeds[mode]))
        cut.append(dict(acoustic_cut_offset_deg=delta, mode_index=mode,
                        **evaluate(changed, rf_mhz=frequency)))
    length = []
    for delta in (-e.dimension_tolerance_mm, e.dimension_tolerance_mm):
        changed = copy.deepcopy(report)
        changed["inputs"]["length_mm"] = s["length_mm"]+delta
        if changed["inputs"]["length_mm"] <= 0:
            length.append(dict(interaction_thickness_offset_mm=delta, evaluation_valid=False,
                               throughput_efficiency_fraction=None, longitudinal_mismatch_m_inv=None))
        else:
            length.append(dict(interaction_thickness_offset_mm=delta, **evaluate(changed)))
    dimension_geometry = []
    for field in ("crystal_length_mm", "crystal_width_mm", "crystal_height_mm"):
        for sign in (-1, 1):
            perturbed = asdict(e)
            perturbed[field] += sign*e.dimension_tolerance_mm
            geometry = crystal_geometry(report, DesignSettings(**perturbed))
            dimension_geometry.append(dict(parameter=field, offset_mm=sign*e.dimension_tolerance_mm,
                                           feasible=geometry["feasible"], failures=geometry["failures"],
                                           minimum_optical_clearance_mm=geometry["minimum_optical_clearance_mm"]))
    return _finite(dict(fixed_rf_power_w=power, fixed_nominal_rf_mhz=frequency,
                        rf_carrier_plus_minus_0p1_percent=rf,
                        internal_alignment=alignment, acoustic_cut=cut,
                        interaction_thickness=length, crystal_dimensions=dimension_geometry,
                        method="Explicit fixed-design perturbations; no RF retuning or reoptimization. Cut perturbations recompute acoustic eigenmode and coupling.",
                        limitations="Independent perturbations, not a worst-case Monte Carlo/yield analysis. L perturbations assume the declared interaction slab is physically realized; the active-area geometry does not establish its field cutoff. Dimension cases reconstruct the nominal facet layout; machining errors at a fixed facet position are not bounded. Internal alignment is not external mount-angle tolerance; transducer resonance, face-cut refraction/coupling, beam diffraction and thermal changes omitted."))


BLOCKERS = [
    "Archive the clarified RF return-loss requirement and qualified acceptance protocol: old 1000:1 versus invalid 1:1000 and provisional VSWR <=1.5 are superseded, not measured match evidence.",
    "Qualify crystal lot, orientation by X-ray/Laue, material tensor, optical activity, absorption and photoelastic dispersion.",
    "Supplier review of actual clipped solid, edge/chamfer allowance, clear apertures, cut/parallelism and surface inspection at operating wavelength.",
    "Precise independent primary-source numerical fabrication-capability verification and written TeO2-specific supplier acceptance missing; proposed allocations are not qualified capabilities.",
    "Design LiNbO3 transducer cut, polarization/mode, thickness/resonance, lateral dimensions and electrode stack; none are determined by bulk TeO2 phase matching.",
    "Specify and qualify the requested vacuum cold-welding metallized bond stack: materials, thickness, acoustic impedance, void fraction, cleanliness, vacuum, pressure, temperature, process and electrode interconnect are unresolved; bare oxide surfaces are not assumed cold-weldable.",
    "Measure complex transducer impedance versus RF, temperature and power; synthesize and qualify matching network, feedthrough and 50-ohm reference plane.",
    "Thermal model and continuous 100% duty test at selected RF are mandatory regardless of requested duty: conversion loss, absorber loss, optical absorption, thermal contact and allowable crystal/bond temperature unresolved.",
    "Establish physical interaction boundaries and effective coupling length: the declared model slab is not automatically realized by the larger acoustic active rectangle; coupling outside it can change the RF power response.",
    "Validate finite optical/acoustic beam overlap, acoustic diffraction/walkoff, absorber placement/reflections and multichannel/multiorder coupling.",
    "Coating prescription and measured polarization-dependent optical transmission/Fresnel loss unresolved.",
    "Independent measured RF detuning, efficiency, extinction, spectral linewidth, rise/fall time, damage and temperature performance required; no values invented.",
    "Independent engineering review, supplier-approved toleranced drawing, assembly process, inspection records and traceable acceptance evidence required before fabrication release.",
]


def _blocker_records(messages, e):
    plans = [
        ("Archive the clarified", "rf-requirement-protocol",
         "Obtain a signed RF requirement; measure calibrated complex S11 with port 2 terminated in 50 ohms across the specified [f_min, f_max] band, at operational power and temperature.",
         f"Archive the clarified replacement: RF return loss >= {e.return_loss_target_db:g} dB, reflected power fraction <=10^(-target_dB/10), and VSWR <=1.065:1. Calibrated results including uncertainty pass; old invalid/reversed VSWR is not acceptance."),
        ("Qualify crystal", "material-lot",
         "Review lot certificate, X-ray/Laue orientation and polarization-dependent optical/acoustic characterization at the requested wavelength.",
         "Accept only traceable lot data and tensor/loss/optical-activity uncertainty bounds approved for this model and operating point; numerical loss and uncertainty limits remain TBD."),
        ("Supplier review", "supplier-drawing",
         "Supplier reviews final directed-face CAD, chamfers, ceramic clearance and clear apertures; inspect dimensions, face normals, interferometric figure and bandwidth-defined roughness.",
         "All supplier-approved dimension, orientation, flatness, surface-quality and aperture requirements pass documented inspection. Provisional requirements and pending roughness limits are not purchase approval."),
        ("Precise independent", "primary-source-verification",
         "Retrieve controlled primary-source revisions and exact material-specific numerical capability statements; obtain written TeO2 supplier acceptance of each proposed requirement.",
         "Auditable source excerpts, revision/date and applicability plus written supplier acceptance exist. Attempted URLs and generic optical capabilities alone do not qualify TeO2."),
        ("Design LiNbO3", "transducer",
         "Model and measure transducer piezoelectric cut/mode, thickness/resonance, total ceramic outline, electrode active area, impedance and acoustic field.",
         "Supplier-approved ceramic/electrode drawing fits the mount including borders; measured resonance, field and conversion meet approved carrier +/-0.1% and power requirements. Thickness/cut/conversion pass limits remain TBD."),
        ("Specify and qualify", "bond",
         "Qualify the requested vacuum cold-welding process and compatible metallization; characterize surface preparation/cleanliness, vacuum, pressure, temperature and acoustic/thermal properties. Inspect thickness/voids, thermal cycling and interconnect strain.",
         "A validated metallized-stack drawing and welding process pass approved thickness, void fraction, loss and thermal/mechanical limits without crystal damage. Layer materials/thicknesses, vacuum/pressure/temperature/duration and numerical acceptance limits remain TBD; no weldability of bare TeO2/LiNbO3 is assumed."),
        ("Measure complex", "rf-match",
         "Calibrate vector RF measurements and forward/reflected power at the feedthrough across RF, temperature and selected power; qualify matching components and resonance.",
         f"Measured RF return loss is >={e.return_loss_target_db:g} dB, reflected fraction <=10^(-target_dB/10), and VSWR <=1.065:1 throughout [f_min, f_max], including uncertainty and without ratings violations; port 2 is terminated in 50 ohms and RF/conversion use a consistent reference plane."),
        ("Thermal model", "continuous-duty-thermal",
         "Calculate and measure crystal, transducer, bond and absorber temperatures and drift under continuous 100% duty at selected RF, regardless of requested duty; include optical absorption and heat-sink contact.",
         "All temperatures, drift and stress stay below approved component/application limits for an approved endurance duration and uncertainty budget. Temperature limits, drift bounds and duration are unresolved; pulsed duty alone cannot qualify release."),
        ("Establish physical", "physical-interaction-length",
         "Map acoustic intensity along the incident/diffracted rays throughout the entire crystal, including outside the declared slab. Qualify field localization/apodization or solve the full spatial coupled-wave problem, and measure the RF-power response.",
         "Approved physical field boundaries or an experimentally qualified integrated coupling model support the selected interaction length and power response within approved error bounds. Active-area fit alone does not establish the slab length or nominal efficiency; error limits and localization design are unresolved."),
        ("Validate finite", "finite-field",
         "Map optical/acoustic fields, walkoff and reflection/diffraction; inspect the as-built geometry including ceramic borders/chamfers and perturbation corners; compare measured diffraction with the conditional model.",
         "Declared bounded Gaussian power, all as-built geometric margins, acoustic coverage and application-approved efficiency/error limits pass. No claim of containing infinite tails or guaranteed modeled overlap."),
        ("Coating prescription", "optical-loss",
         "Qualify coating prescription and measure polarization-resolved transmission, absorption and scatter over the required wavelengths/angles and thermal conditions.",
         "Approved transmission, wavefront, scatter and damage limits pass measured acceptance; assumed unity transmission and preliminary scratch/dig alone are not evidence. Numerical limits remain TBD."),
        ("Independent measured", "application-performance",
         "Measure efficiency, RF detuning, extinction, spectral linewidth, modulation rise/fall and damage/temperature behavior with calibrated instruments.",
         "User-approved numerical performance limits and measurement uncertainty pass at the required operating point; linewidth, extinction, timing and temperature requirements remain unresolved."),
        ("Independent engineering", "release-review",
         "Independent engineering/QA review checks supplier-approved drawings, assembly process, all qualification records and every blocker acceptance criterion.",
         "Every blocker has traceable passing evidence and signed independent review. No status field, approval toggle or simulated result constitutes fabrication release."),
        ("Geometry infeasible", "geometry",
         "Reconstruct and inspect a physically fitted solid and finite optical/acoustic envelopes including approved machining allowances.",
         "All declared finite-envelope geometry constraints are nonnegative and supplier-approved clear apertures, ceramic/bond margins and chamfers fit as built."),
        ("Selected RF power", "rf-ceiling",
         "Check calibrated selected forward source power at the agreed reference plane, including uncertainty and power-control overshoot.",
         f"Selected RF power and uncertainty/overshoot remain within the approved engineering ceiling {e.rf_ceiling_w:g} W, strictly below 1 W."),
    ]
    records = []
    for message in messages:
        match = next((item for item in plans if message.startswith(item[0])), None)
        if match is None:
            raise ValueError("Missing verification plan for release blocker")
        _, identifier, verification, criterion = match
        records.append(dict(id=identifier, statement=message, status="open",
                            verification=verification, acceptance_criterion=criterion,
                            evidence=None))
    return records


def _optical_apertures(geometry):
    result = []
    for footprint in geometry["optical_footprints"]:
        name = footprint["name"]
        if "entrance" not in name and "exit" not in name:
            continue
        center = np.array(footprint["center_local_mm"])
        extent = np.array(footprint["envelope_half_extents_local_mm"])
        result.append(dict(name=name, datum="A" if "entrance" in name else "B",
                           center_local_mm=center.tolist(),
                           ellipse_axis1_local_mm=footprint["envelope_axis1_local_mm"],
                           ellipse_axis2_local_mm=footprint["envelope_axis2_local_mm"],
                           y_interval_mm=[float(center[1]-extent[1]), float(center[1]+extent[1])],
                           z_interval_mm=[float(center[2]-extent[2]), float(center[2]+extent[2])]))
    return result


def _svg_crystal(geometry, tolerances):
    vertices = np.array(geometry["vertex_coordinates_local_mm"])
    dimensions = geometry["nominal_dimensions_mm"]
    rotation = np.array(geometry["local_to_crystal_4x4"])[:3, :3]
    datums = {"entrance": "A", "exit": "B", "acoustic_mount": "C",
              "side_minus_y": "D", "side_plus_y": "E", "side_minus_z": "F",
              "side_plus_z": "G"}
    apertures = _optical_apertures(geometry)
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="1050" height="510" viewBox="0 0 1050 510">',
           '<rect width="1050" height="510" fill="white"/>',
           '<defs><marker id="direction" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0 0 L6 3 L0 6 Z" fill="black"/></marker></defs>',
           '<text x="20" y="24" font-size="17">TeO2 actual clipped solid — orthographic projections (mm)</text>',
           f'<text x="20" y="47" fill="red">Revision {DESIGN_REVISION} — NOT FOR FABRICATION; chamfers/ceramic borders not included</text>']
    scale = min(230/max(dimensions), 35)
    for panel, (a, b, title) in enumerate(((0, 1, "x-y / top"), (0, 2, "x-z / side"),
                                          (1, 2, "y-z / end"))):
        cx, cy = 175+panel*340, 245

        def pos(p):
            return cx+p[a]*scale, cy-p[b]*scale

        for face in geometry["faces"]:
            polygon = vertices[face["vertex_indices"]]
            coordinates = " ".join(f"{x:.5f},{y:.5f}" for x, y in map(pos, polygon))
            fill = "#f4d6b5" if face["name"] == "acoustic_mount" else "#c6dbea"
            out.append(f'<polygon points="{coordinates}" fill="{fill}" fill-opacity=".25" stroke="#26465c"/>')
            center = np.array(face["centroid_local_mm"])
            normal = np.array(face["outward_normal_local"])
            x, y = pos(center)
            nx, ny = pos(center+normal*0.5)
            datum = datums[face["name"]]
            out.append(f'<path d="M{x:.5f},{y:.5f} L{nx:.5f},{ny:.5f}" stroke="#7a3a00" fill="none" marker-end="url(#direction)" data-datum="{datum}"/>')
            out.append(f'<text x="{x+4:.5f}" y="{y+11:.5f}" font-size="9">{datum}: {escape(face["name"])}</text>')
        for aperture in apertures:
            center = np.array(aperture["center_local_mm"])
            v1 = np.array(aperture["ellipse_axis1_local_mm"])
            v2 = np.array(aperture["ellipse_axis2_local_mm"])
            points = [center+v1*math.cos(t)+v2*math.sin(t)
                      for t in np.linspace(0, 2*math.pi, 65)]
            coordinates = " ".join(f"{x:.5f},{y:.5f}" for x, y in map(pos, points))
            out.append(f'<polyline points="{coordinates}" stroke="#a71925" stroke-width="1.3" fill="none"/>')
        for axis_index, (axis, label, color) in enumerate(zip(
                np.eye(3), ("[100]", "[010]", "[001]"), ("#a71925", "#167b2c", "#204ba0"))):
            projected = rotation.T@axis
            x0, y0 = cx-100, 88
            x1, y1 = x0+35*projected[a], y0-35*projected[b]
            out.append(f'<path d="M{x0:.5f},{y0:.5f} L{x1:.5f},{y1:.5f}" stroke="{color}" marker-end="url(#direction)" data-crystal-axis="{label}" data-view="{"xyz"[a]}{"xyz"[b]}"/>')
            out.append(f'<text x="{x1+2:.5f}" y="{y1-2-11*axis_index:.5f}" font-size="10" fill="{color}">{label}</text>')
        for i, vertex in enumerate(vertices):
            x, y = pos(vertex)
            out.append(f'<text x="{x+3:.4f}" y="{y-3:.4f}" font-size="9">V{i}</text>')
        x0, x1 = cx-dimensions[a]*scale/2, cx+dimensions[a]*scale/2
        yy = cy+dimensions[b]*scale/2+25
        xx = x0-16
        y0, y1 = cy-dimensions[b]*scale/2, cy+dimensions[b]*scale/2
        out.extend([f'<path d="M{x0:.5f},{yy-7:.5f} V{yy+7:.5f} M{x0:.5f},{yy:.5f} H{x1:.5f} M{x1:.5f},{yy-7:.5f} V{yy+7:.5f}" stroke="black" fill="none"/>',
                    f'<text x="{cx:.5f}" y="{yy+23:.5f}" text-anchor="middle" font-size="13">{dimensions[a]:.6g} mm nominal {"xyz"[a]}</text>',
                    f'<path d="M{xx-5:.5f},{y0:.5f} H{xx+5:.5f} M{xx:.5f},{y0:.5f} V{y1:.5f} M{xx-5:.5f},{y1:.5f} H{xx+5:.5f}" stroke="black" fill="none"/>',
                    f'<text x="{xx-8:.5f}" y="{cy:.5f}" transform="rotate(-90 {xx-8:.5f} {cy:.5f})" text-anchor="middle" font-size="11">{dimensions[b]:.6g} mm nominal {"xyz"[b]}</text>',
                    f'<text x="{cx:.5f}" y="385" text-anchor="middle" font-size="13">{escape(title)}; nominal {"xyz"[b]}={dimensions[b]:.6g} mm</text>'])
    out.append('<text x="20" y="422" font-size="12">x = incident k; z = interaction-plane normal. Crystal axes use computed transform projections.</text>')
    out.append(f'<text x="20" y="444" font-size="12">Actual volume {geometry["volume_mm3"]:.6g} mm^3; scale {scale:.6g} SVG units/mm; orange = acoustic facet.</text>')
    out.append('<text x="20" y="468" font-size="12">Red curves: bounded optical apertures. Brown lines: projected outward normals. Orange: acoustic facet C.</text>')

    def vector(value):
        return "("+", ".join(f"{x:.7g}" for x in value)+")"

    mount = geometry["acoustic_mount"]
    slab = geometry["interaction_slab"]
    angles = geometry["cut_angles"]
    operating = geometry["operating_point"]
    rows = [
        f"REVISION {DESIGN_REVISION}; design_status=preliminary_concept; no supplier/independent approval evidence.",
        f"Operating point: vacuum wavelength={operating['wavelength_nm']:.7g} nm; RF carrier={operating['rf_mhz']:.7g} MHz; optical shift sign={operating['diffraction_frequency_shift_sign']:+d}; acoustic eigenmode index={operating['acoustic_mode_index']}.",
        f"RF interface requirement: return loss strictly >{operating['return_loss_target_db']:.7g} dB. Matching network and calibrated acceptance evidence absent.",
        "DATUMS: A entrance (outward -incident k), B parallel exit (+incident k), C acoustic mount (outward -q).",
        "Origin: nominal blank center. Local x,y,z and crystal X=[100],Y=[010],Z=[001] are right-handed.",
        "FACE PLANES: outward unit n; material satisfies n dot local_position <= offset_mm. Polygons use vertex IDs below.",
    ]
    for face in geometry["faces"]:
        rows.extend([
            f"{datums[face['name']]} {face['name']}: local n={vector(face['outward_normal_local'])}; offset={face['plane_offset_mm']:.7g} mm; area={face['area_mm2']:.7g} mm2",
            f"    crystal n={vector(face['outward_normal_crystal'])}; polygon={face['vertex_indices']}",
        ])
    rows += ["TRANSFORM local -> crystal (rotation rows; translation=0 mm):"]
    rows += [vector(row) for row in rotation]
    rows += [
        "CUT ANGLES (degrees; directed vectors, not merely acute cut):",
        f"  Basal azimuth={angles['azimuth_deg']:.7g}; q to +[001]={angles['acoustic_normal_to_plus_c_deg']:.7g}; acoustic plane to C acute={angles['acoustic_face_plane_to_c_acute_deg']:.7g}",
        f"  Incident signed theta={angles['incident_wavevector_signed_to_c_deg']:.7g}; diffracted signed theta={angles['diffracted_wavevector_signed_to_c_deg']:.7g}",
        f"  Polarizations i={geometry['optical_polarizations']['incident']}, d={geometry['optical_polarizations']['diffracted']}; q crystal={vector(geometry['acoustic_wave_normal_unit_crystal'])}",
        f"INTERACTION SLAB: x={slab['x_interval_mm']} mm; thickness={slab['thickness_mm']:.7g} mm, versus crystal face separation={slab['crystal_optical_face_separation_mm']:.7g} mm.",
        "  Slab is a model assumption: physical RF-field cutoff/effective coupling length unresolved; active area can couple outside it.",
        f"  Bounded optical prism lower={vector(slab['bounded_optical_envelope_box_local_mm'][0])}, upper={vector(slab['bounded_optical_envelope_box_local_mm'][1])} mm",
        f"  Optical-envelope volume upper bound={slab['bounded_optical_volume_upper_bound_mm3']:.7g} mm3; NOT Gaussian intensity or overlap volume.",
        f"TRANSDUCER ACTIVE AREA: {mount['wavefront_dimensions_mm']} mm, {mount['active_acoustic_area_mm2']:.7g} mm2; center={vector(mount['center_local_mm'])} mm.",
        f"  Active width axis={vector(mount['width_axis_local'])}; height axis={vector(mount['height_axis_local'])}; inward normal={vector(mount['inward_normal_local'])}.",
        "  TOTAL CERAMIC footprint/thickness, electrode border, bond/interconnect and edge clearance: UNRESOLVED; active rectangle is not ceramic outline.",
        f"BOUNDED APERTURES: Gaussian 1/e2 diameter={geometry['beam']['diameter_mm']:.7g} mm; radius={geometry['beam']['bounded_radius_mm']:.7g} mm; tail power={geometry['beam']['excluded_power_fraction']:.7g}.",
        "  Ellipse = center + axis1*cos(t) + axis2*sin(t). Infinite Gaussian tails are excluded, not contained.",
    ]
    for aperture in apertures:
        rows += [
            f"  {aperture['datum']} {aperture['name']}: center={vector(aperture['center_local_mm'])}; axis1={vector(aperture['ellipse_axis1_local_mm'])}; axis2={vector(aperture['ellipse_axis2_local_mm'])} mm",
            f"    y range={vector(aperture['y_interval_mm'])}; z range={vector(aperture['z_interval_mm'])} mm",
        ]
    rows += [
        "PROVISIONAL TOLERANCES — engineering allocations only; source retrieval UNVERIFIED, supplier acceptance missing:",
        f"  Blank dimensions +/-{tolerances['dimension']['value']:.7g} mm; acoustic crystallographic cut +/-{tolerances['acoustic_cut']['value']:.7g} deg.",
        f"  Substrate crystal-frame orientation +/-{tolerances['substrate_orientation']['value']:.7g} deg; absolute optical-face normal +/-{tolerances['optical_face_orientation']['value']:.7g} deg.",
        "  Entrance/exit parallelism <=30 arcsec; scratch/dig proposal 20-10; chamfer proposal 0.1 mm NOT included in geometry.",
        "  Surface figure proposal lambda/10 PV SURFACE HEIGHT at 632.8 nm (63.28 nm), over declared bounded optical aperture union.",
        "  Roughness metric: areal RMS height Sq. Numerical roughness limit, sampling/filter bandwidth, scan size/instrument: PENDING.",
        "  Flatness aperture guard band/chamfer allowance, total ceramic fit and as-built tolerance stack require a revised supplier drawing.",
        "  Acceptance: inspect each approved dimension/normal/aperture; qualify RF match, continuous 100% duty thermal and measured optical performance.",
        "  Per-value rationale, attempted primary-source URLs and individual blocker verification/acceptance plans are in crystal_spec.json.",
        "VERTICES (mm), sufficient with face polygons and rotation above to reconstruct the convex solid:",
    ]
    for i, point in enumerate(vertices):
        rows.append(f"V{i}: local={vector(point)}; crystal={vector(rotation@point)}")
    y = 500
    for row in rows:
        for line in textwrap.wrap(row, width=135) or [""]:
            out.append(f'<text x="20" y="{y}" font-size="12" font-family="monospace">{escape(line)}</text>')
            y += 18
    height = y+25
    out[0] = f'<svg xmlns="http://www.w3.org/2000/svg" width="1050" height="{height}" viewBox="0 0 1050 {height}">'
    out[1] = f'<rect width="1050" height="{height}" fill="white"/>'
    out.append("</svg>")
    return "\n".join(out)+"\n"


def _svg_assembly():
    return """<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="420" viewBox="0 0 1000 420">
<rect width="1000" height="420" fill="white"/>
<text x="20" y="30" font-size="20">ASSEMBLY SCHEMATIC — NOT TO SCALE — NOT FOR FABRICATION</text>
<rect x="395" y="125" width="190" height="160" fill="#d9e9f3" stroke="black"/>
<text x="420" y="205">Clipped TeO2 crystal</text>
<path d="M50 195 H390 M590 180 L925 150 M590 220 H925" stroke="#a22" fill="none"/>
<text x="35" y="175">Collimated, polarized input</text>
<text x="700" y="135">First order / analyzer / stop</text>
<text x="700" y="245">Zero order / beam dump</text>
<rect x="420" y="85" width="140" height="34" fill="#f4d6b5" stroke="black"/>
<text x="365" y="72">LiNbO3 / electrodes / bond — TBD</text>
<path d="M100 340 H345 V100 H415" stroke="#267" fill="none"/>
<text x="40" y="365">RF source / gating / matching / feedthrough — TBD</text>
<rect x="410" y="290" width="160" height="30" fill="#ddd" stroke="black"/>
<text x="590" y="310">Absorber / heat sink / compliant mount — TBD</text>
<text x="20" y="400">Connections are functional only; physical positions follow actual facet geometry after supplier and thermal review.</text>
</svg>
"""


def write_design_package(report, engineering, prefix):
    """Write six strict, auditable concept artifacts; return their file paths."""
    e = _engineering(engineering)
    geometry = crystal_geometry(report, e)
    r, inputs = report["optimum"], report["inputs"]
    power = float(r["rf_selected_w"])
    gamma = None if e.measured_vswr == 0 else (e.measured_vswr-1)/(e.measured_vswr+1)
    accepted = None if gamma is None else power*(1-gamma**2)
    reflection_limit = 10**(-e.return_loss_target_db/10)
    gamma_limit = math.sqrt(reflection_limit)
    vswr_limit = None if gamma_limit >= 1 else (1+gamma_limit)/(1-gamma_limit)
    measured_return_loss = None if gamma is None or gamma == 0 else -20*math.log10(gamma)
    meets_return_loss = (None if gamma is None else e.measured_vswr <= vswr_limit)
    blockers = BLOCKERS.copy()
    if not geometry["feasible"]:
        blockers.insert(0, "Geometry infeasible: "+"; ".join(geometry["failures"]))
    if power > e.rf_ceiling_w+1e-12:
        blockers.insert(0, "Selected RF power exceeds engineering RF ceiling")
    reference = dict(
        dimension=dict(value=e.dimension_tolerance_mm, unit="mm", meaning="bilateral ±, nominal blank dimensions",
                       selection_rationale="Small dimensional allocation relative to the 10 mm concept blank; explicit +/- perturbations check sensitivity. No process capability or manufacturing yield is established.",
                       source_status="UNVERIFIED engineering-proposed requirement, not supplier capability",
                       attempted_primary_sources=["Newport fabrication capabilities matrix", "Edmund optical specification notation"],
                       evidence="Proposed engineering allocation; requires supplier capability confirmation, not a measured TeO2 tolerance"),
        acoustic_cut=dict(value=e.cut_tolerance_deg, unit="deg", meaning="± crystallographic acoustic wave-normal tilt",
                          selection_rationale="Trial crystallographic-cut allocation evaluated by recalculating acoustic eigenmode and fixed-RF diffraction. The resulting sensitivity, not the number of printed digits, must inform a later approved tolerance.",
                          source_status="UNVERIFIED engineering-proposed requirement, not supplier capability",
                          attempted_primary_sources=["Crys-Teh supplier qualification", "Newport fabrication capabilities matrix"],
                          evidence="Proposed allocation tested by explicit acoustic-cut perturbations; not a supplier specification"),
        internal_alignment=dict(value=e.alignment_tolerance_deg, unit="deg", meaning="± internal optical wavevector, not external mount angle",
                                selection_rationale="Fine internal-angle perturbation for an initial adjustment budget; response is calculated but external mount/refraction sensitivity and an acceptable efficiency loss are not yet approved.",
                                source_status="UNVERIFIED engineering-proposed requirement, not demonstrated alignment capability",
                                attempted_primary_sources=["Edmund optical specification notation"],
                                evidence="Proposed allocation evaluated by fixed-design response; no experimental acceptance evidence"),
        surface_quality=dict(proposed_scratch_dig="20-10", proposed_flatness="lambda/10 at 632.8 nm",
                             proposed_parallelism_arcsec=30, edge_chamfer_mm=0.1,
                             per_value_proposals=dict(
                                 scratch_dig=dict(value="20-10", unit="scratch/dig designation, not two defect widths",
                                                  selection_rationale="Initial inspection target for supplier quotation. Scattering/extinction budgets have not established its adequacy, and no TeO2 supplier capability is asserted.",
                                                  source_status="UNVERIFIED engineering-proposed requirement",
                                                  attempted_primary_sources=["Edmund optical specification notation", "Crys-Teh supplier qualification"]),
                                 surface_flatness=dict(value=0.1, unit="waves PV surface height over approved clear aperture",
                                                       reference_wavelength_nm=632.8,
                                                       equivalent_surface_height_pv_nm=63.28,
                                                       selection_rationale="Initial polished-surface figure allocation at a stated metrology wavelength. It is NOT a transmitted-wavefront or reflected-wavefront guarantee; optical-path sensitivity of both faces remains to be qualified.",
                                                       source_status="UNVERIFIED engineering-proposed requirement",
                                                       attempted_primary_sources=["Newport fabrication capabilities matrix", "Edmund optical specification notation"]),
                                 parallelism=dict(value=30, unit="arcsec maximum entrance/exit normal deviation",
                                                  selection_rationale="Trial small angular allocation, about 0.00833 degrees, comparable in scale to the proposed internal alignment budget. These are not interchangeable tolerances; face-refraction sensitivity and supplier feasibility remain unresolved.",
                                                  source_status="UNVERIFIED engineering-proposed requirement",
                                                  attempted_primary_sources=["Newport fabrication capabilities matrix", "Crys-Teh supplier qualification"]),
                                 edge_chamfer=dict(value=0.1, unit="mm nominal edge break",
                                                  selection_rationale="Trial edge-break allowance intended to reduce handling chips. Not included in the current solid or aperture checks; geometry and transducer clearances must be revised before use.",
                                                  source_status="UNVERIFIED engineering-proposed requirement",
                                                  attempted_primary_sources=["Crys-Teh supplier qualification", "Edmund optical specification notation"])),
                             status="Illustrative procurement proposals only; supplier qualification, physical aperture allowance and wavefront sensitivity unresolved",
                             source="No verified TeO2-specific capability source. Do not treat these values as certified or Crys-teh claims."),
        definitions_source=dict(title="MIL-PRF-13830B, Optical Components for Fire Control Instruments; General Specification Governing the Manufacture, Assembly, and Inspection of",
                                url="https://quicksearch.dla.mil/",
                                use="Scratch/dig terminology only; does not establish TeO2 supplier capability or justify the proposed numerical tolerances",
                                verification="Exact controlled revision and applicable clauses must be retrieved and agreed with supplier before purchase."),
        primary_source_verification=dict(status="UNRESOLVED — official source retrieval failed DNS; no numerical supplier capability is asserted",
                                         sources_to_verify=[
                                             dict(title="Crys-Teh supplier qualification",
                                                  url="https://crys-teh.com/",
                                                  use="Obtain material-, geometry- and lot-specific written cut/polish/bond capability and quotation",
                                                  numerical_claim_verified=False),
                                             dict(title="Newport Custom Optics Fabrication Capabilities Matrix",
                                                  url="https://www.newport.com/medias/sys_master/images/images/hd8/h1a/8989202415646/Newport-Custom-Optics-Fabrication-Capabilities-Matrix.pdf",
                                                  use="Generic fabrication capability comparison only; material-specific qualification still needed",
                                                  numerical_claim_verified=False),
                                             dict(title="Edmund Optics, Understanding Optical Specifications",
                                                  url="https://www.edmundoptics.com/knowledge-center/application-notes/optics/understanding-optical-specifications/",
                                                  use="Specification definitions and notation; not a TeO2 capability guarantee",
                                                  numerical_claim_verified=False)]))
    for key, meaning, rationale in (
            ("substrate_orientation", "± absolute crystallographic frame registration per approved rotational datum",
             "Trial registration allocation at the same scale as the acoustic cut budget. X-ray/Laue metrology and coupled three-dimensional optical/acoustic sensitivity must establish a final tolerance; this is not tested by the single tilt perturbation alone."),
            ("optical_face_orientation", "± absolute entrance/exit face normal relative to the prescribed crystal frame; common absolute error separate from parallelism",
             "Trial optical-face cut allocation. Common face-normal error can change refraction and alignment; the 30 arcsec relative parallelism proposal is separately tighter. No supplier capability or full face-cut sensitivity qualification is implied.")):
        reference[key] = dict(value=e.cut_tolerance_deg, unit="deg", meaning=meaning,
                              selection_rationale=rationale,
                              source_status="UNVERIFIED engineering-proposed requirement",
                              attempted_primary_sources=["Crys-Teh supplier qualification", "Newport fabrication capabilities matrix"],
                              evidence=None)
    surface = reference["surface_quality"]["per_value_proposals"]
    reference["substrate_orientation"]["angular_error_definition"] = (
        "Maximum geodesic rotation angle between nominal and measured right-handed crystal frames; "
        "not independently summed Euler-angle tolerances.")
    reference["optical_face_orientation"]["angular_error_definition"] = (
        "Angle between the prescribed and measured directed unit face normals <= proposed value; "
        "+/- denotes either tilt sign, not separately compounded pitch/yaw limits.")
    surface["parallelism"]["angular_error_definition"] = (
        "acos(-n_A dot n_B), since entrance/exit outward normals are antiparallel, <= proposed arcsec.")
    surface["surface_flatness"]["reference_aperture"] = dict(
        frame="local drawing frame, mm", datum_faces=["A", "B"],
        region="Union of explicitly bounded optical ellipse footprints on each face",
        footprints=_optical_apertures(geometry), additional_guard_band_mm=None,
        edge_exclusion_mm=None,
        qualification="Guard band, edge exclusion, chamfers, interferometer fitting/removal convention and full required clear aperture must be supplier-approved; surface figure is not transmitted wavefront.")
    surface["roughness"] = dict(
        value=None, unit="nm areal RMS height Sq",
        metric="Areal RMS surface height Sq; roughness separated from figure by an explicitly approved spatial filtering/sampling bandwidth",
        spatial_bandwidth_mm_inv=None, sampling_pitch_um=None, scan_size_mm=None,
        instrument=None, filter_and_form_removal=None,
        selection_rationale="A bandwidth-defined RMS height metric is required to distinguish roughness/scatter from long-scale surface figure. No scattering model, retrieved capability or user acceptance budget justifies inventing a numerical limit or bandwidth.",
        source_status="UNVERIFIED requirement definition; numerical limit and measurement bandwidth PENDING",
        attempted_primary_sources=["Edmund optical specification notation", "Crys-Teh supplier qualification"],
        evidence=None)
    blocker_records = _blocker_records(blockers, e)
    bom = [
        dict(item="TeO2 crystal", quantity=1, specification="Actual clipped solid and directed crystal frame in crystal_spec.json", status="supplier qualification required"),
        dict(item="LiNbO3 transducer", quantity=1, specification="Cut/thickness/electrodes/resonance TBD; acoustic active area from geometry", status="release blocker"),
        dict(item="Vacuum cold-welded acoustic bond and interconnect", quantity=1, specification="Compatible metallization/layers, thickness, surface preparation, vacuum/pressure/temperature and void/loss/damage tests TBD; no bare-oxide weldability assumed", status="release blocker"),
        dict(item="RF matching / feedthrough", quantity=1, specification=f"Carrier {r['rf_mhz']:.6g} MHz, selected source power {power:.6g} W; measured impedance needed", status="release blocker"),
        dict(item="Acoustic absorber / thermal mount", quantity=1, specification="Location, loss, thermal resistance and compliance TBD", status="release blocker"),
        dict(item="AR coating / optical analyzer / beam stops", quantity=1, specification=f"{inputs['wavelength_nm']:.6g} nm, polarized inputs/outputs; transmission qualification TBD", status="release blocker"),
    ]
    acceptance = [
        dict(test="Dimensional and crystallographic inspection", criterion="Supplier-approved dimensions, directed face normals, cut and optical parallelism; finite-envelope aperture recheck including chamfers", evidence=None),
        dict(test="Optical surface and transmission", criterion="Agree scratch/dig, surface-PV flatness at 632.8 nm over the declared aperture plus approved guard band, bandwidth-defined areal Sq roughness, coating and polarization; measure transmission. Roughness limits/bandwidth pending.", evidence=None),
        dict(test="RF impedance / return loss", criterion=f"Measure calibrated input S11 with port 2 terminated in 50 Ω across [f_min, f_max] (numerical endpoints to be agreed), power and temperature. Return loss must be >={e.return_loss_target_db:g} dB, reflected forward-power fraction <={reflection_limit:g}, and VSWR <=1.065:1, including uncertainty; legacy VSWR <=1.5 proposal is superseded.", evidence=None),
        dict(test="Conditional diffraction model", criterion="Compare measured efficiency with modeled nominal and ±0.1% RF at fixed selected power; measurement uncertainty and pass threshold must be approved", evidence=None),
        dict(test="Thermal / endurance", criterion=f"Run continuous 100% duty regardless of requested {e.duty_cycle_percent:g}% duty; measure crystal/bond/transducer/absorber temperatures and drift. Approved safe maxima, uncertainty and endurance duration TBD.", evidence=None),
        dict(test="Dynamic, spectral and extinction", criterion="User-approved linewidth, rise/fall and extinction requirements TBD; no inferred guarantees", evidence=None),
    ]
    spec = dict(status="PRELIMINARY / CONCEPT — NOT FOR FABRICATION",
                design_status="preliminary_concept", allowed_design_states=DESIGN_STATES,
                design_revision=DESIGN_REVISION,
                design_status_reasons=[b["statement"] for b in blocker_records],
                geometry=geometry,
                proposed_tolerances=reference,
                datums=dict(A="entrance optical plane, outward -incident k", B="exit plane parallel A",
                            C="acoustic mounting plane, inward +q", origin="nominal blank center",
                            crystal_axes="X=[100], Y=[010], Z=[001]; all vectors directed"),
                manufacturing_sequence=["Qualify blank and crystallographic frame", "Supplier-approved cut and mount facet",
                                        "Recompute apertures after approved chamfers and tolerances",
                                        "Polish and inspect optical faces", "Qualify coating", "Bond qualified transducer and inspect",
                                        "Assemble absorber/thermal mount", "Measure RF impedance, tune matching and perform acceptance"],
                release_blockers=blocker_records)
    paths = {name: str(prefix)+suffix for name, suffix in
             (("design_json", "_design.json"), ("design_report", "_design.md"),
              ("crystal_spec_json", "_crystal_spec.json"), ("crystal_spec_text", "_crystal_spec.txt"),
              ("crystal_drawing", "_crystal.svg"), ("assembly_schematic", "_assembly.svg"))}
    manifest = dict(schema_version="1.1", status=spec["status"], design_status="preliminary_concept",
                    allowed_design_states=DESIGN_STATES, design_revision=DESIGN_REVISION,
                    design_status_reasons=spec["design_status_reasons"],
                    fabrication_release=False,
                    release_policy="Evidence-based only; no approval toggle. Independent review, qualified drawings and measured acceptance records are missing.",
                    engineering_inputs=asdict(e), optimizer_inputs=inputs, material=report["material"],
                    optimizer_version=report.get("version"), optimizer_model_warnings=report.get("warnings", []),
                    literature_audit=report.get("literature_audit", "See AUDIT.md and VALIDATION.md; no experimental device qualification."),
                    selected_design=r, geometry=geometry, modeled_performance=_performance(report, e),
                    rf=dict(selected_forward_source_power_w=power, ceiling_w=e.rf_ceiling_w,
                            ceiling_met=power <= e.rf_ceiling_w+1e-12,
                            duty_cycle_percent=e.duty_cycle_percent,
                            average_forward_source_power_w=power*e.duty_cycle_percent/100,
                            target_vswr=e.vswr_target, measured_vswr=None if gamma is None else e.measured_vswr,
                            target_vswr_status="Superseded provisional setting, not the current RF acceptance criterion",
                            return_loss_target_db=e.return_loss_target_db,
                            return_loss_comparison="greater than or equal to",
                            reflected_power_fraction_limit=reflection_limit,
                            reflected_power_fraction_comparison="less than or equal to",
                            s11_db_upper_limit=-e.return_loss_target_db,
                            equivalent_vswr_limit=vswr_limit,
                            stated_vswr_limit=1.065,
                            stated_vswr_comparison="less than or equal to",
                            rf_match_frequency_band="[f_min, f_max]; numerical endpoints to be agreed",
                            port2_termination_ohm=50,
                            return_loss_requirement_confirmed=True,
                            return_loss_requirement_basis="Clarified replacement RF requirement recorded in README; a specified requirement is not measured matching evidence",
                            measured_reflection_magnitude=gamma, net_accepted_rf_w=accepted,
                            measured_return_loss_db=measured_return_loss,
                            supplied_vswr_is_zero_reflection_limit=gamma == 0 if gamma is not None else None,
                            return_loss_meets_target_from_supplied_vswr=meets_return_loss,
                            measured_vswr_meets_stated_limit=None if gamma is None else e.measured_vswr <= 1.065,
                            measured_point_meets_rf_limits=None if gamma is None else
                                meets_return_loss and e.measured_vswr <= 1.065,
                            reflected_power_limit_at_selected_rf_w=power*reflection_limit,
                            match_evidence="User-supplied single-point VSWR only; calibrated band-sweep evidence with port 2 terminated in 50 ohms is still required" if gamma is not None else "No measured match supplied; calibrated band sweep with port 2 terminated in 50 ohms is required",
                            measured_vswr_is_band_compliance_evidence=False,
                            target_reflected_fraction=reflection_limit,
                            legacy_vswr_reflected_fraction=((e.vswr_target-1)/(e.vswr_target+1))**2,
                            net_accepted_note="Paccepted=Pforward*(1-|Gamma|^2) only at the same calibrated reference plane. Target never substituted for measurement; cable loss unprovided.",
                            requirement_history=dict(original="1000:1", subsequent="1:1000",
                                                     standard="VSWR >= 1", reversal_valid=False,
                                                     resolved=True, superseded=True,
                                                     resolution="Current requirement is inclusive RF return loss >= the configured dB threshold across [f_min, f_max], with port 2 terminated in 50 ohms")),
                    conditional_assumptions=dict(rf_to_acoustic_conversion=inputs["conversion"],
                                                 optical_transmission=inputs["optical_transmission"],
                                                 overlap_intensity_factor=inputs["overlap"],
                                                 interaction_slab_physically_realized=None,
                                                 interaction_length_status="Declared model thickness; active-area fit does not prove field localization or full-crystal coupled-wave response",
                                                 efficiency_time_basis="Steady RF-on two-wave response at selected forward power, not duty-cycle-averaged optical efficiency. Envelope transients and repetition timing are unmodeled.",
                                                 interpretation="User model assumptions, not measured efficiency/transmission/match. Conversion denominator must be calibrated; accepted RF is not silently reapplied to avoid double counting."),
                    unmodeled_performance=dict(spectral_linewidth_nm=None, rise_time_ns=None,
                                               extinction_ratio=None, allowable_temperature_c=None,
                                               thermal_rise_c=None, impedance_ohm=None),
                    proposed_tolerances=reference, bom=bom, acceptance=acceptance,
                    release_blockers=blocker_records, evidence=dict(measured_acceptance_records=[],
                                                           supplier_drawing_approval=None,
                                                           independent_review=None),
                    artifacts=paths)
    nominal = manifest["modeled_performance"]["rf_carrier_plus_minus_0p1_percent"][1]
    lines = ["# TeO2 engineering concept package", "", "**PRELIMINARY / CONCEPT — NOT FOR FABRICATION**", "",
             "## Selected operating point",
             f"- Wavelength: {inputs['wavelength_nm']:.6g} nm; RF: {r['rf_mhz']:.6g} MHz.",
             f"- Selected fixed RF source power: {power:.6g} W; ceiling: {e.rf_ceiling_w:.6g} W; duty: {e.duty_cycle_percent:g}%.",
             f"- Conditional nominal throughput: {nominal['throughput_efficiency_fraction']}; geometry feasible: {geometry['feasible']}.",
             "- Modeled conversion, transmission and overlap are conditional assumptions, not measured evidence.",
             f"- Measured VSWR: {manifest['rf']['measured_vswr']}; net accepted RF: {accepted} W (null/None means unavailable).",
             f"- RF input-match requirement: return loss >={e.return_loss_target_db:g} dB (S11 <= {-e.return_loss_target_db:g} dB), reflected fraction <={reflection_limit:g}, and VSWR <=1.065:1 over [f_min, f_max], with port 2 terminated in 50 Ω.",
             f"- Exact VSWR limit corresponding to {e.return_loss_target_db:g} dB return loss is {vswr_limit:.9g}:1; the stated 1.065:1 cap is slightly stricter. Numerical band endpoints remain to be agreed.",
             "- Standard VSWR is >=1. Original 1000:1 then invalid reversal 1:1000 and provisional <=1.5:1 are superseded by the clarified return-loss requirement, not measured matching evidence.",
             "", "## Geometry and beam envelope",
             f"- Nominal blank: {e.crystal_length_mm:g} x {e.crystal_width_mm:g} x {e.crystal_height_mm:g} mm; clipped volume {geometry['volume_mm3']:.6g} mm^3.",
             "- Entrance and parallel exit are normal to incident k; extraordinary rays are computed separately from k.",
             f"- 1/e^2 diameter {e.beam_diameter_mm:g} mm; bounded radius {geometry['beam']['bounded_radius_mm']:.6g} mm; excluded Gaussian power {e.gaussian_tail_fraction:g}.",
             "- Gaussian tails are infinite; finite-envelope containment is not a full finite-beam overlap integral.",
             "- Air exit conserves tangential dimensionless k exactly and uses sqrt(1-|b_t|^2), with explicit TIR flags, not scalar extraordinary Snell.",
             "- Numeric vertices, face polygons/outward normals, halfspaces, transforms, beam footprints and acoustic group-walkoff checks are in the JSON spec.",
             f"- Interaction slab: {geometry['interaction_slab']['thickness_mm']:g} mm, not the {e.crystal_length_mm:g} mm crystal length. Bounded optical-prism volume upper bound: {geometry['interaction_slab']['bounded_optical_volume_upper_bound_mm3']:.6g} mm^3, not an intensity/overlap volume.",
             "- The slab is a model assumption, not a verified RF-field cutoff. The active area can cause coupling outside it; localization or a full spatial coupled-wave calculation and measured RF-power response are release blockers.",
             f"- Fitted acoustic active area: {inputs['acoustic_width_mm']:g} x {inputs['acoustic_height_mm']:g} mm. Total ceramic, electrode border and bond/interconnect footprint remain unresolved.",
             "", "## Fixed-design performance and tolerances",
             "- RF carrier +/-0.1% is evaluated at the same chosen power, wavelength and cut; no retuning.",
             "- Acoustic-cut perturbations recompute the acoustic eigenmode and coupling; alignment and interaction thickness are perturbed explicitly.",
             "- All proposed numerical tolerances are allocations, not demonstrated TeO2 manufacturing capability. Surface proposals require supplier agreement.",
             "- No spectral linewidth, measured rise time, extinction, temperature limit or transducer resonance is invented.",
             "```json", json.dumps(_finite(manifest["modeled_performance"]), indent=2, allow_nan=False), "```",
             "", "## Proposed tolerances (not released)",
             "```json", json.dumps(reference, indent=2, allow_nan=False), "```",
             "", "## Bill of materials"]
    lines += [f"- {item['item']} ({item['quantity']}): {item['specification']} — {item['status']}." for item in bom]
    lines += ["", "## Acceptance plan (all evidence missing)"]
    lines += [f"- {item['test']}: {item['criterion']}." for item in acceptance]
    lines += ["", "## Release blockers"]
    for item in blocker_records:
        lines += [f"- **{item['id']}**: {item['statement']}",
                  f"  - Verification: {item['verification']}",
                  f"  - Acceptance: {item['acceptance_criterion']}; evidence: missing."]
    lines += ["", "No approval toggle constitutes release. Independent review and traceable supplier/measurement evidence are mandatory."]
    text = ["PRELIMINARY / CONCEPT — NOT FOR FABRICATION", f"Revision {DESIGN_REVISION}; design_status=preliminary_concept",
            "Coordinates and dimensions in mm.",
            "Nominal blank dimensions: "+", ".join(f"{v:.9g}" for v in geometry["nominal_dimensions_mm"]),
            "Final solid: convex halfspace intersection; acoustic facet removes blank material.",
            "Vertices in drawing frame (x along incident k):"]
    text += [f"V{i}: "+", ".join(f"{v:.9g}" for v in point) for i, point in
             enumerate(geometry["vertex_coordinates_local_mm"])]
    text += ["Faces (outward local normal, plane n.x=offset, polygon vertex IDs):"]
    text += [f"{face['name']}: n={face['outward_normal_local']}; offset={face['plane_offset_mm']:.9g}; vertices={face['vertex_indices']}"
             for face in geometry["faces"]]
    text += ["Local-to-crystal homogeneous transform:"]+[str(row) for row in geometry["local_to_crystal_4x4"]]
    text += ["Interaction slab:", json.dumps(geometry["interaction_slab"], indent=2, allow_nan=False),
             "Acoustic active area and ceramic status:",
             json.dumps(geometry["acoustic_mount"], indent=2, allow_nan=False)]
    text += ["Proposed tolerances (not supplier-qualified):", json.dumps(reference, indent=2, allow_nan=False),
             "Release blockers and verification/acceptance plans:",
             json.dumps(blocker_records, indent=2, allow_nan=False)]
    payloads = dict(design_json=json.dumps(_finite(manifest), indent=2, allow_nan=False)+"\n",
                    design_report="\n".join(lines)+"\n",
                    crystal_spec_json=json.dumps(_finite(spec), indent=2, allow_nan=False)+"\n",
                    crystal_spec_text="\n".join(text)+"\n",
                    crystal_drawing=_svg_crystal(geometry, reference), assembly_schematic=_svg_assembly())
    for name, path in paths.items():
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(payloads[name], encoding="utf-8")
    return paths
