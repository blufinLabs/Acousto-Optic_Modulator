#!/usr/bin/env python3
"""Preliminary bulk TeO2 acousto-optic cut optimizer (Python 3.10+).

Uses the two-wave framework and parallel-tangent design discussed in Chang's
1976 review, with explicit numerical extensions. See README.md and AUDIT.md.
Units inside the model: SI; command-line units: nm, mm, MHz, W, degrees.
The optimization is in one specified crystal plane, not a global 3-D search.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar

VERSION = "0.2.0"
SCOPE = "preliminary non-gyrotropic plane-wave bulk Bragg model"


@dataclass(frozen=True)
class Material:
    """Room-temperature constants; Mys et al. (2019), section 2.

    Elastic coefficients in GPa; photoelastic coefficients dimensionless,
    referenced to 632.8 nm. Their spectral dispersion is NOT modeled.
    Optical indices use our two-oscillator refit of Uchida (1971), Table I;
    see fit_dispersion.py for the primary data and reproducible fit.
    """
    rho: float = 5990.0
    c11: float = 53.2
    c12: float = 48.6
    c13: float = 21.2
    c33: float = 108.5
    c44: float = 24.4
    c66: float = 55.2
    p11: float = 0.0074
    p12: float = 0.187
    p13: float = 0.340
    p31: float = 0.0905
    p33: float = 0.240
    p44: float = -0.170
    p66: float = -0.0463

    def validate(self):
        if not all(math.isfinite(v) for v in asdict(self).values()):
            raise ValueError("Material constants must be finite")
        if self.rho <= 0:
            raise ValueError("Density must be positive")
        c = self.stiffness()
        if np.min(np.linalg.eigvalsh(c)) <= 0:
            raise ValueError("Stiffness matrix must be positive definite")

    def stiffness(self):
        a, b, c, d, e, f = (getattr(self, k) * 1e9 for k in
                           ("c11", "c12", "c13", "c33", "c44", "c66"))
        return np.array([[a, b, c, 0, 0, 0], [b, a, c, 0, 0, 0],
                         [c, c, d, 0, 0, 0], [0, 0, 0, e, 0, 0],
                         [0, 0, 0, 0, e, 0], [0, 0, 0, 0, 0, f]])

    def photoelastic(self):
        return np.array([[self.p11, self.p12, self.p13, 0, 0, 0],
                         [self.p12, self.p11, self.p13, 0, 0, 0],
                         [self.p31, self.p31, self.p33, 0, 0, 0],
                         [0, 0, 0, self.p44, 0, 0],
                         [0, 0, 0, 0, self.p44, 0],
                         [0, 0, 0, 0, 0, self.p66]])


def indices(wavelength_nm: float, allow_extrapolation=False):
    """Two-term Sellmeier; extrapolation beyond 1000 nm requires opt-in."""
    upper = 4500 if allow_extrapolation else 1000
    if not math.isfinite(wavelength_nm) or not 400 <= wavelength_nm <= upper:
        raise ValueError("Use 400-1000 nm, or explicitly enable dispersion extrapolation up to 4500 nm")
    w2 = (wavelength_nm / 1000) ** 2
    resonance1, resonance2 = 1.239841984/9.24, 1.239841984/4.70
    no = math.sqrt(1 + 2.584013829506934*w2/(w2-resonance1**2) +
                   1.1563319953188333*w2/(w2-resonance2**2))
    ne = math.sqrt(1 + 2.8361380477666955*w2/(w2-resonance1**2) +
                   1.528349297324733*w2/(w2-resonance2**2))
    return no, ne


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def angle(a, b):
    return float(np.degrees(np.arccos(np.clip(np.dot(unit(a), unit(b)), -1, 1))))


def christoffel(q, m: Material):
    x, y, z = q
    a, b, c, d, e, f = (getattr(m, k)*1e9 for k in
                         ("c11", "c12", "c13", "c33", "c44", "c66"))
    return np.array([[a*x*x+f*y*y+e*z*z, (b+f)*x*y, (c+e)*x*z],
                     [(b+f)*x*y, f*x*x+a*y*y+e*z*z, (c+e)*y*z],
                     [(c+e)*x*z, (c+e)*y*z, e*(x*x+y*y)+d*z*z]])


def acoustic_modes(q, m):
    """Increasing phase speed; use displacement vectors, not assumed shear."""
    eigenvalues, displacements = np.linalg.eigh(christoffel(unit(q), m))
    return np.sqrt(eigenvalues/m.rho), displacements


def acoustic_group(q, displacement, speed, m):
    """Gradient of homogeneous omega(K); tracks polarization at crossings."""
    h = 1e-6
    tangent_gradient = np.empty(3)
    for j in range(3):
        dq = np.eye(3)[j]*h
        vp, up = acoustic_modes(q+dq, m)
        vm, um = acoustic_modes(q-dq, m)
        ip = np.argmax(np.abs(up.T @ displacement))
        im = np.argmax(np.abs(um.T @ displacement))
        tangent_gradient[j] = (vp[ip]-vm[im])/(2*h)
    tangent_gradient -= np.dot(tangent_gradient, q)*q
    return speed*q + tangent_gradient


def impermeability(q, displacement, m):
    """Delta(1/n^2) per k_acoustic * displacement amplitude.

    Engineering strain [Sxx,Syy,Szz,2Syz,2Sxz,2Sxy]; output has
    [Bxx,Byy,Bzz,Byz,Bxz,Bxy], without doubling optical off-diagonals.
    """
    s = (np.outer(q, displacement)+np.outer(displacement, q))/2
    sv = np.array([s[0, 0], s[1, 1], s[2, 2],
                   2*s[1, 2], 2*s[0, 2], 2*s[0, 1]])
    b = m.photoelastic() @ sv
    return np.array([[b[0], b[5], b[4]], [b[5], b[1], b[3]],
                     [b[4], b[3], b[2]]])


def rotational_impermeability(q, displacement, no, ne):
    """Kinematic crystal rotation: delta B = W B0 - B0 W.

    W is the antisymmetric part of the displacement gradient. This explicit
    extension is derived from R B0 R.T; it is not a formula transcribed from
    Chang. A calibrated effective coefficient may already contain this effect.
    """
    h = np.outer(displacement, q)
    w = (h-h.T)/2
    b0 = np.diag([1/no**2, 1/no**2, 1/ne**2])
    return w@b0-b0@w


def flux_projection(direction, d, no, ne, normal):
    """Normalized Poynting flux through fixed interaction boundary.

    With |D|=1, E=B0 D, H proportional to n*s cross E.
    Return J = n [(s.normal)*|E|² - (s.E)*(E.normal)].
    For isotropic n this reduces to cos(theta)/n³.
    """
    directions, ds = np.atleast_2d(direction), np.atleast_2d(d)
    normals = np.broadcast_to(normal, directions.shape)
    b0 = np.diag([1/no**2, 1/no**2, 1/ne**2])
    es = ds @ b0
    ns = 1/np.sqrt(np.einsum("ni,ij,nj->n", ds, b0, ds))
    return ns*(np.sum(directions*normals, axis=1)*np.sum(es*es, axis=1)-
               np.sum(directions*es, axis=1)*np.sum(es*normals, axis=1))


def two_wave_efficiency(kappa2, mismatch, length, transmission=1.0):
    """Chang Eq.13 functional form with g=kappa² L², delta=Delta k L/2."""
    kappa2, mismatch = np.broadcast_arrays(kappa2, mismatch)
    denominator = kappa2+(mismatch/2)**2
    ratio = np.divide(kappa2, denominator, out=np.zeros_like(denominator, dtype=float),
                      where=denominator > 0)
    return transmission*ratio*np.sin(length*np.sqrt(denominator))**2


@dataclass(frozen=True)
class Settings:
    wavelength_nm: float
    design_mode: str = "efficiency"
    strain_only: bool = False
    allow_dispersion_extrapolation: bool = False
    azimuth_deg: float = 45.0
    tilt_min_deg: float = -30.0
    tilt_max_deg: float = 30.0
    theta_min_deg: float = 10.0
    theta_max_deg: float = 85.0
    rf_min_mhz: float = 20.0
    rf_max_mhz: float = 200.0
    length_mm: float = 5.0
    acoustic_width_mm: float = 5.0
    acoustic_height_mm: float = 1.0
    rf_budget_w: float = 1.0
    conversion: float = 0.5
    overlap: float = 1.0
    optical_transmission: float = 1.0
    loss_db_per_mm_100mhz: float = 0.0
    acoustic_path_mm: float = 0.0
    max_walkoff_deg: float = 85.0
    grid_deg: float = 2.0
    interaction: str = "anisotropic"
    incident: str = "both"
    acoustic_mode: str = "all"

    def validate(self):
        for k, v in asdict(self).items():
            if isinstance(v, (float, int)) and not math.isfinite(v):
                raise ValueError(f"{k} must be finite")
        indices(self.wavelength_nm, self.allow_dispersion_extrapolation)
        if not -80 <= self.tilt_min_deg < self.tilt_max_deg <= 80:
            raise ValueError("Tilt bounds must increase and lie within [-80,80]")
        if not 5 <= self.theta_min_deg < self.theta_max_deg <= 85:
            raise ValueError("Optical angle bounds must lie within [5,85] degrees")
        if not 0 < self.rf_min_mhz <= self.rf_max_mhz:
            raise ValueError("RF bounds must be positive and increasing")
        for k in ("length_mm", "acoustic_width_mm", "acoustic_height_mm", "grid_deg"):
            if getattr(self, k) <= 0:
                raise ValueError(f"{k} must be positive")
        if not 0.25 <= self.grid_deg <= 10:
            raise ValueError("grid_deg must lie within [0.25,10]")
        for k in ("rf_budget_w", "loss_db_per_mm_100mhz", "acoustic_path_mm"):
            if getattr(self, k) < 0:
                raise ValueError(f"{k} cannot be negative")
        for k in ("conversion", "overlap", "optical_transmission"):
            if not 0 < getattr(self, k) <= 1:
                raise ValueError(f"{k} must lie within (0,1]")
        if not 0 < self.max_walkoff_deg <= 85:
            raise ValueError("max_walkoff_deg must lie within (0,85]")
        if self.interaction not in ("anisotropic", "isotropic", "all"):
            raise ValueError("Unknown interaction type")
        if self.incident not in ("o", "e", "both"):
            raise ValueError("Unknown incident polarization")
        if self.acoustic_mode not in ("slow", "middle", "fast", "all"):
            raise ValueError("Unknown acoustic mode")
        if self.design_mode not in ("efficiency", "wide-aperture"):
            raise ValueError("Unknown design mode")
        if self.design_mode == "wide-aperture" and self.interaction != "anisotropic":
            raise ValueError("wide-aperture mode requires anisotropic interaction")


class Model:
    def __init__(self, settings: Settings, material: Material | None = None):
        self.s = settings
        self.m = material or Material()
        self.s.validate()
        self.m.validate()
        self.no, self.ne = indices(settings.wavelength_nm, settings.allow_dispersion_extrapolation)
        phi = np.radians(settings.azimuth_deg)
        self.t = np.array([np.cos(phi), np.sin(phi), 0.0])
        self.z = np.array([0.0, 0.0, 1.0])
        self.o = np.cross(self.z, self.t)
        # Dimensionless wave vector b = k * lambda/(2*pi).
        self.surfaces = {"o": np.eye(3)/self.no**2,
                         "e": np.diag([1/self.ne**2, 1/self.ne**2, 1/self.no**2])}

    def direction(self, theta_deg):
        a = np.radians(np.atleast_1d(theta_deg))
        return np.sin(a)[:, None]*self.t + np.cos(a)[:, None]*self.z

    def optical(self, directions, pol):
        directions = np.atleast_2d(directions)
        if pol == "o":
            return np.full(len(directions), self.no), np.tile(self.o, (len(directions), 1))
        n = 1/np.sqrt(np.einsum("ni,ij,nj->n", directions,
                                self.surfaces["e"], directions))
        # D is perpendicular to k; E is generally not for an extraordinary wave.
        d = unit(np.cross(np.tile(self.o, (len(directions), 1)), directions))
        return n, d

    def branch(self, alpha_deg, theta_deg, mode, pi, pd, root):
        """Vectorized exact geometric phase matching, then approximate coupling.

        root=-1/+1 selects a quadratic root; its sign does NOT prescribe order.
        Zero-frequency roots, backward waves and near-axis eigenstates excluded.
        """
        alpha = math.radians(float(alpha_deg))
        q = np.cos(alpha)*self.t + np.sin(alpha)*self.z
        speeds, displacements = acoustic_modes(q, self.m)
        v, u = speeds[mode], displacements[:, mode]
        vg = acoustic_group(q, u, v, self.m)
        walkoff = angle(q, vg)
        si = self.direction(theta_deg)
        ni, di = self.optical(si, pi)
        bi = ni[:, None]*si
        ad = self.surfaces[pd]
        a = q @ ad @ q
        b = 2*(bi @ ad @ q)
        c = np.einsum("ni,ij,nj->n", bi, ad, bi)-1
        discr = b*b-4*a*c
        valid = discr >= -1e-13
        # Stable quadratic roots avoid cancellation for the near-forward branch.
        sq = np.sqrt(np.maximum(discr, 0))
        far = -0.5*(b+np.copysign(sq, b)) / a
        near = np.divide(c, a*far, out=np.zeros_like(c), where=np.abs(far) > 1e-15)
        low, high = np.minimum(near, far), np.maximum(near, far)
        Q = low if root == -1 else high
        bd = bi + Q[:, None]*q
        nd = np.linalg.norm(bd, axis=1)
        sd = unit(bd)
        _, dd = self.optical(sd, pd)
        theta_d = np.degrees(np.arctan2(sd @ self.t, sd[:, 2]))
        freq = np.abs(Q)*v/(self.s.wavelength_nm*1e-9)/1e6
        strain_b = impermeability(q, u, self.m)
        rotation_b = rotational_impermeability(q, u, self.no, self.ne)
        p_strain = np.einsum("ni,ij,nj->n", dd, strain_b, di)
        p_rotation = np.einsum("ni,ij,nj->n", dd, rotation_b, di)
        p = p_strain if self.s.strain_only else p_strain+p_rotation
        m2 = ni**3 * nd**3 * p**2 / (self.m.rho*v**3)
        ri = unit(bi @ self.surfaces[pi])
        rd = unit(bd @ ad)
        ray_separation = np.degrees(np.arccos(np.clip(np.sum(ri*rd, axis=1), -1, 1)))
        parallel_cross = np.cross(ri, rd) @ self.o
        ji = flux_projection(si, di, self.no, self.ne, si)
        jd = flux_projection(sd, dd, self.no, self.ne, si)
        flux_factor = np.divide(1, ni**3*nd**3*ji*jd,
                                out=np.zeros_like(ni), where=(ji > 0) & (jd > 0))
        loss_db = (self.s.loss_db_per_mm_100mhz * self.s.acoustic_path_mm *
                   (freq/100)**2)
        delivery = self.s.conversion*self.s.overlap*10**(-loss_db/10)
        L = self.s.length_mm*1e-3
        area = self.s.acoustic_width_mm*self.s.acoustic_height_mm*1e-6
        gain = np.pi**2*L**2*m2*flux_factor*delivery/(2*(self.s.wavelength_nm*1e-9)**2*area)
        # Geometric validity is independent of configurable RF limits.
        valid &= (np.abs(Q) > 1e-10) & (sd[:, 2] > 0) & (np.sum(si*sd, axis=1) > 0)
        valid &= ((np.abs(theta_d) >= self.s.theta_min_deg) &
                  (np.abs(theta_d) <= self.s.theta_max_deg))
        valid &= walkoff <= self.s.max_walkoff_deg+1e-7
        valid &= np.isfinite(gain) & (gain > 1e-14)
        return dict(valid=valid, gain=gain, rf_mhz=freq, m2=m2, p=p,
                    p_strain=p_strain, p_rotation=p_rotation, flux_factor=flux_factor,
                    parallel_cross=parallel_cross,
                    theta_d=theta_d, ray_separation=ray_separation, v=float(v),
                    q=q, u=u, vg=vg, walkoff=walkoff, ni=ni, nd=nd,
                    Q=Q, loss_db=loss_db, delivery=delivery,
                    residual=np.abs(np.einsum("ni,ij,nj->n", bd, ad, bd)-1))

    def record(self, alpha, theta, mode, pi, pd, root):
        r = self.branch(alpha, [theta], mode, pi, pd, root)
        f = float(r["rf_mhz"][0])
        if not r["valid"][0] or not self.s.rf_min_mhz-1e-5 <= f <= self.s.rf_max_mhz+1e-5:
            return None
        if self.s.design_mode == "wide-aperture" and abs(r["parallel_cross"][0]) > 1e-7:
            return None
        g = float(r["gain"][0])
        peak = math.pi**2/(4*g)
        chosen = min(self.s.rf_budget_w, peak)
        ideal = math.sin(math.sqrt(g*chosen))**2
        fixed = math.sin(math.sqrt(g*self.s.rf_budget_w))**2
        return dict(acoustic_tilt_from_basal_deg=float(alpha),
                    acoustic_normal_to_c_axis_deg=90-float(alpha),
                    acoustic_face_plane_to_c_axis_deg=abs(float(alpha)),
                    incident_wavevector_signed_to_c_deg=float(theta),
                    incident_wavevector_to_c_axis_deg=abs(float(theta)),
                    diffracted_wavevector_signed_to_c_deg=float(r["theta_d"][0]),
                    diffracted_wavevector_to_c_axis_deg=abs(float(r["theta_d"][0])),
                    azimuth_deg=self.s.azimuth_deg, incident_polarization=pi,
                    diffracted_polarization=pd, acoustic_mode_index=int(mode),
                    acoustic_phase_speed_m_s=r["v"], acoustic_walkoff_deg=r["walkoff"],
                    acoustic_normal_xyz=r["q"].tolist(), acoustic_displacement_xyz=r["u"].tolist(),
                    acoustic_group_velocity_xyz_m_s=r["vg"].tolist(), rf_mhz=f,
                    frequency_shift_sign=int(np.sign(r["Q"][0])),
                    ni=float(r["ni"][0]), nd=float(r["nd"][0]),
                    effective_photoelastic=float(r["p"][0]),
                    strain_photoelastic_contribution=float(r["p_strain"][0]),
                    rotational_contribution=float(r["p_rotation"][0]),
                    optical_flux_correction=float(r["flux_factor"][0]),
                    parallel_tangent_residual=float(r["parallel_cross"][0]),
                    klein_cook_q_diagnostic=float(2*np.pi*self.s.wavelength_nm*1e-9*self.s.length_mm*1e-3 /
                        (r["ni"][0]*(r["v"]/(f*1e6))**2)),
                    m2_s3_kg=float(r["m2"][0]), gain_per_rf_w=g,
                    rf_first_peak_w=peak, rf_selected_w=chosen,
                    intrinsic_efficiency=ideal,
                    throughput_efficiency=ideal*self.s.optical_transmission,
                    efficiency_at_exact_budget=fixed,
                    acoustic_power_at_interaction_w=chosen*float(r["delivery"][0]),
                    acoustic_loss_db=float(r["loss_db"][0]),
                    optical_ray_separation_deg=float(r["ray_separation"][0]),
                    dimensionless_surface_residual=float(r["residual"][0]),
                    acoustic_transit_diagnostic_ns=1e6*self.s.acoustic_height_mm/r["v"],
                    quadratic_root=int(root))

    def wide_geometry(self, ordinary_theta_deg, mode, pi):
        """Exact circle/ellipse parallel tangents; returns a phase-matched branch.

        tan(theta_e)=(ne/no)^2*tan(theta_o); this follows from equality of
        the optical surface normals, rather than small-birefringence expansion.
        """
        theta_o = float(ordinary_theta_deg)
        theta_e = math.degrees(math.atan((self.ne/self.no)**2*math.tan(math.radians(theta_o))))
        so, se = self.direction([theta_o])[0], self.direction([theta_e])[0]
        neff = float(self.optical([se], "e")[0][0])
        delta = neff*se-self.no*so
        if np.dot(delta, self.t) < 0:
            delta = -delta
        q = unit(delta)
        alpha = math.degrees(math.atan2(float(q[2]), float(q@self.t)))
        theta, theta_d, pd = (theta_o, theta_e, "e") if pi == "o" else (theta_e, theta_o, "o")
        roots = [(root, self.branch(alpha, [theta], mode, pi, pd, root)) for root in (-1, 1)]
        root, r = min(roots, key=lambda x: abs(float(x[1]["theta_d"][0])-theta_d))
        return alpha, theta, pd, root, r


def response(report, *, wavelength_nm=None, rf_mhz=None, alignment_offset_deg=None, rf_power_w=None):
    """Fixed-cut response with transverse momentum conserved.

    Solve the diffracted index surface along the nominal incident propagation
    coordinate, and use the longitudinal mismatch in the detuned two-wave
    coupled-mode solution. The interaction boundary normal stays fixed.
    Alignment offsets are INTERNAL optical wavevector angles, in the chosen
    plane; neither the acoustic cut nor its frequency is retuned automatically.
    """
    s = Settings(**report["inputs"])
    m = Material(**report["material"])
    optimum = report["optimum"]
    wavelength_nm = s.wavelength_nm if wavelength_nm is None else wavelength_nm
    rf_mhz = optimum["rf_mhz"] if rf_mhz is None else rf_mhz
    alignment_offset_deg = 0.0 if alignment_offset_deg is None else alignment_offset_deg
    rf_power_w = optimum["rf_selected_w"] if rf_power_w is None else rf_power_w
    ws, fs, offsets, powers = np.broadcast_arrays(wavelength_nm, rf_mhz, alignment_offset_deg, rf_power_w)
    eta, mismatches = np.zeros(ws.shape), np.full(ws.shape, np.nan)
    q = np.asarray(optimum["acoustic_normal_xyz"])
    u = np.asarray(optimum["acoustic_displacement_xyz"])
    v = optimum["acoustic_phase_speed_m_s"]
    phi = math.radians(s.azimuth_deg)
    t = np.array([math.cos(phi), math.sin(phi), 0])
    o = np.cross([0, 0, 1], t)
    strain_b = impermeability(q, u, m)
    reference_theta = math.radians(optimum["incident_wavevector_signed_to_c_deg"])
    normal = math.sin(reference_theta)*t+math.cos(reference_theta)*np.array([0, 0, 1])
    L = s.length_mm*1e-3
    area = s.acoustic_width_mm*s.acoustic_height_mm*1e-6
    pi, pd = optimum["incident_polarization"], optimum["diffracted_polarization"]
    for idx in np.ndindex(ws.shape):
        w, f, power = float(ws[idx]), float(fs[idx]), float(powers[idx])
        if not all(math.isfinite(x) for x in (w, f, power, float(offsets[idx]))) or f <= 0 or power < 0:
            raise ValueError("Response inputs must be finite, with positive RF frequency and nonnegative power")
        no, ne = indices(w, s.allow_dispersion_extrapolation)
        surfaces = {"o": np.eye(3)/no**2, "e": np.diag([1/ne**2, 1/ne**2, 1/no**2])}
        theta = math.radians(optimum["incident_wavevector_signed_to_c_deg"]+float(offsets[idx]))
        si = math.sin(theta)*t+math.cos(theta)*np.array([0, 0, 1])
        # Refuse extrapolation of the non-gyrotropic eigenstates toward the optic axis.
        if not 5 <= abs(math.degrees(theta)) <= 85:
            eta[idx] = np.nan
            continue
        ni = 1/math.sqrt(si @ surfaces[pi] @ si)
        di = o if pi == "o" else unit(np.cross(o, si))
        db = strain_b if s.strain_only else strain_b+rotational_impermeability(q, u, no, ne)
        target = ni*si + optimum["frequency_shift_sign"]*w*1e-9*f*1e6/v*q
        beta_target = float(target @ normal)
        transverse = target-beta_target*normal
        ad = surfaces[pd]
        a, b, c = normal @ ad @ normal, 2*normal @ ad @ transverse, transverse @ ad @ transverse-1
        discriminant = b*b-4*a*c
        if discriminant < 0:
            continue
        roots = [(-b+math.sqrt(discriminant))/(2*a), (-b-math.sqrt(discriminant))/(2*a)]
        forward_roots = [r for r in roots if r > 0]
        if not forward_roots:
            continue
        beta = min(forward_roots, key=lambda x: abs(x-beta_target))
        bd = transverse+beta*normal
        nd = np.linalg.norm(bd)
        sd = bd/nd
        theta_d = abs(math.degrees(math.atan2(sd @ t, sd[2])))
        if not 5 <= theta_d <= 85:
            eta[idx] = np.nan
            continue
        dd = o if pd == "o" else unit(np.cross(o, sd))
        peff = float(dd @ db @ di)
        ji = float(flux_projection(si, di, no, ne, normal)[0])
        jd = float(flux_projection(sd, dd, no, ne, normal)[0])
        if ji <= 0 or jd <= 0:
            continue
        coupling_m2 = peff**2/(m.rho*v**3*ji*jd)
        loss_db = s.loss_db_per_mm_100mhz*s.acoustic_path_mm*(f/100)**2
        intensity = power*s.conversion*s.overlap*10**(-loss_db/10)/area
        kappa2 = math.pi**2/(w*1e-9)**2*coupling_m2*intensity/2
        mismatch = 2*math.pi/(w*1e-9)*(beta_target-beta)
        mismatches[idx] = mismatch
        eta[idx] = two_wave_efficiency(kappa2, mismatch, L, s.optical_transmission)
    return eta, mismatches


def response_sweeps(report, samples=1601):
    """Four illustrative fixed-design sweeps; returns arrays and actual settings."""
    s = Settings(**report["inputs"])
    r = report["optimum"]
    w = s.wavelength_nm
    wl_lo, wl_hi = max(400, w-w*0.015), min(4500 if s.allow_dispersion_extrapolation else 1000, w+w*0.015)
    f = r["rf_mhz"]
    axes = {
        "wavelength_nm": np.linspace(wl_lo, wl_hi, samples),
        "rf_mhz": np.linspace(max(0.001, f-f*0.1), f+f*0.1, samples),
        "alignment_offset_deg": np.linspace(-1, 1, samples),
        "rf_power_w": np.linspace(0, max(2.5*r["rf_first_peak_w"], s.rf_budget_w*1.1), samples),
    }
    result = {}
    for key, x in axes.items():
        y, mismatch = response(report, **{key: x})
        result[key] = dict(x=x, efficiency=y, mismatch_m_inv=mismatch)
    return result


def write_response_graphs(report, prefix):
    import matplotlib.pyplot as plt
    sweeps = response_sweeps(report)
    labels = ["Vacuum wavelength (nm)", "RF drive frequency (MHz)",
              "Internal optical alignment offset (degrees)", "RF input power (W)"]
    titles = ["Wavelength sensitivity — fixed RF", "RF-frequency sensitivity — fixed wavelength",
              "Angular sensitivity — fixed wavelength and RF", "RF-power response — matched operating point"]
    r = report["optimum"]
    markers = [report["inputs"]["wavelength_nm"], r["rf_mhz"], 0, r["rf_selected_w"]]
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (key, sweep), xlabel, title, marker in zip(axs.flat, sweeps.items(), labels, titles, markers):
        x, y = sweep["x"], sweep["efficiency"]
        ax.plot(x, y*100, color="#176da0", lw=1.5)
        ax.axvline(marker, color="#a34822", ls="--", lw=1, label="selected operating point")
        ax.set(xlabel=xlabel, ylabel="Diffraction efficiency / incident power (%)", title=title, ylim=(0, 105))
        ax.grid(alpha=0.25)
        if key == "rf_power_w":
            ax.axvline(r["rf_first_peak_w"], color="black", ls=":", label="first theoretical peak")
            ax.legend(fontsize=8)
        # CSV contains actual data, not a screenshot of a plot.
        with Path(str(prefix)+"_"+key+".csv").open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([key, "throughput_efficiency_fraction", "longitudinal_mismatch_m_inv"])
            writer.writerows(zip(x, y, sweep["mismatch_m_inv"]))
    fig.suptitle(f"TeO₂ {report['inputs']['wavelength_nm']:g} nm — {report['inputs']['design_mode']} design (v{VERSION})", fontsize=14)
    fig.text(0.5, 0.015,
             f"Fixed acoustic cut; nominal RF power {r['rf_selected_w']:.4g} W except power sweep. "
             "Fixed boundary; internal angles; no retuning. Optical activity omitted.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(str(prefix)+"_performance.png", dpi=180)
    fig.savefig(str(prefix)+"_performance.pdf")
    plt.close(fig)


def sampled_angles(lo, hi, step):
    return np.linspace(lo, hi, int(math.ceil((hi-lo)/step))+1)


def optimize(settings: Settings, material: Material | None = None):
    """Coarse deterministic scan + constrained multistart refinement.

    Maximize first-lobe coupling per delivered RF watt. With tunable RF power
    below a budget this also maximizes reachable efficiency; once saturated,
    choose the cut requiring least RF power. No certified global optimum.
    """
    model = Model(settings, material)
    if settings.design_mode == "wide-aperture":
        candidates, scan = optimize_wide(model)
        return build_report(settings, model, candidates, scan)
    modes = [0, 1, 2] if settings.acoustic_mode == "all" else [
        {"slow": 0, "middle": 1, "fast": 2}[settings.acoustic_mode]]
    pols = [(a, b) for a in "oe" for b in "oe"
            if (settings.incident == "both" or a == settings.incident)
            and (settings.interaction == "all" or
                 (a != b) == (settings.interaction == "anisotropic"))]
    positive = sampled_angles(settings.theta_min_deg, settings.theta_max_deg, settings.grid_deg)
    thetas = np.concatenate((-positive[::-1], positive))
    alphas = sampled_angles(settings.tilt_min_deg, settings.tilt_max_deg, settings.grid_deg)
    candidates = []
    scan = []
    fixed_frequency = settings.rf_min_mhz == settings.rf_max_mhz
    for mode in modes:
        for pi, pd in pols:
            for root in (-1, 1):
                for side in (-1, 1):
                    ts = thetas[thetas*side > 0]
                    seeds = []
                    for alpha in alphas:
                        r = model.branch(alpha, ts, mode, pi, pd, root)
                        admissible = r["valid"] & (r["rf_mhz"] >= settings.rf_min_mhz) & (
                            r["rf_mhz"] <= settings.rf_max_mhz)
                        # Fixed RF: seed near the equality constraint, then solve it below.
                        seed_valid = r["valid"] if fixed_frequency else admissible
                        score = r["gain"] / (1+((r["rf_mhz"]-settings.rf_min_mhz)/10)**2) if fixed_frequency else r["gain"]
                        for i in np.flatnonzero(seed_valid):
                            seeds.append((float(score[i]), float(alpha), float(ts[i])))
                        if fixed_frequency:
                            differences = r["rf_mhz"]-settings.rf_min_mhz
                            for j in range(len(ts)-1):
                                if differences[j]*differences[j+1] > 0:
                                    continue
                                try:
                                    th = brentq(lambda x: float(model.branch(alpha, [x], mode, pi, pd, root)["rf_mhz"][0])-
                                                settings.rf_min_mhz, ts[j], ts[j+1], xtol=1e-10)
                                except ValueError:
                                    continue
                                rec = model.record(alpha, th, mode, pi, pd, root)
                                if rec:
                                    candidates.append(rec)
                                    seeds.append((rec["gain_per_rf_w"], float(alpha), th))
                                    scan.append(dict(tilt_deg=float(alpha), mode=mode, interaction=pi+pd,
                                                     gain_per_rf_w=rec["gain_per_rf_w"],
                                                     m2_s3_kg=rec["m2_s3_kg"], rf_mhz=rec["rf_mhz"]))
                        for i in np.flatnonzero(admissible):
                            rec = model.record(alpha, ts[i], mode, pi, pd, root)
                            if rec:
                                candidates.append(rec)
                        if np.any(admissible):
                            i = np.argmax(np.where(admissible, r["gain"], -np.inf))
                            scan.append(dict(tilt_deg=float(alpha), mode=mode, interaction=pi+pd,
                                             gain_per_rf_w=float(r["gain"][i]),
                                             m2_s3_kg=float(r["m2"][i]), rf_mhz=float(r["rf_mhz"][i])))
                    seeds.sort(reverse=True)
                    # Spatially distinct starts reduce sensitivity to one local basin.
                    chosen_seeds = []
                    for seed in seeds:
                        if all(np.hypot(seed[1]-s[1], seed[2]-s[2]) > 4 for s in chosen_seeds):
                            chosen_seeds.append(seed)
                        if len(chosen_seeds) == 6:
                            break
                    bounds = [(settings.tilt_min_deg, settings.tilt_max_deg),
                              (settings.theta_min_deg, settings.theta_max_deg) if side == 1 else
                              (-settings.theta_max_deg, -settings.theta_min_deg)]

                    def data(x):
                        return model.branch(x[0], [x[1]], mode, pi, pd, root)

                    def objective(x):
                        rr = data(x)
                        return -float(rr["gain"][0])

                    def geometry(x):
                        rr = data(x)
                        td = abs(float(rr["theta_d"][0]))
                        return np.array([td-settings.theta_min_deg, settings.theta_max_deg-td,
                                         settings.max_walkoff_deg-rr["walkoff"]])

                    def rf_constraints(x):
                        f = float(data(x)["rf_mhz"][0])
                        return np.array([f-settings.rf_min_mhz, settings.rf_max_mhz-f])

                    constraints = [{"type": "ineq", "fun": geometry}]
                    if fixed_frequency:
                        constraints.append({"type": "eq", "fun": lambda x:
                                            float(data(x)["rf_mhz"][0])-settings.rf_min_mhz})
                    else:
                        constraints.append({"type": "ineq", "fun": rf_constraints})
                    for _, alpha, theta in chosen_seeds:
                        result = minimize(objective, [alpha, theta], method="SLSQP", bounds=bounds,
                                          constraints=constraints,
                                          options={"maxiter": 150, "ftol": 1e-10})
                        if result.success:
                            rec = model.record(*result.x, mode, pi, pd, root)
                            if rec:
                                candidates.append(rec)
    return build_report(settings, model, candidates, scan)


def optimize_wide(model):
    """Optimize on the exact parallel-tangent geometry family (one parameter)."""
    s = model.s
    modes = [0, 1, 2] if s.acoustic_mode == "all" else [{"slow": 0, "middle": 1, "fast": 2}[s.acoustic_mode]]
    pols = "oe" if s.incident == "both" else s.incident
    candidates, scan = [], []
    fixed = s.rf_min_mhz == s.rf_max_mhz
    for mode in modes:
        for pi in pols:
            for side in (-1, 1):
                lo, hi = (s.theta_min_deg, s.theta_max_deg) if side == 1 else (-s.theta_max_deg, -s.theta_min_deg)
                angles = sampled_angles(lo, hi, s.grid_deg)

                def raw(theta_o):
                    return model.wide_geometry(theta_o, mode, pi)

                def record(theta_o):
                    alpha, theta, pd, root, _ = raw(theta_o)
                    if not s.tilt_min_deg-1e-8 <= alpha <= s.tilt_max_deg+1e-8:
                        return None
                    if not s.theta_min_deg-1e-8 <= abs(theta) <= s.theta_max_deg+1e-8:
                        return None
                    return model.record(alpha, theta, mode, pi, pd, root)

                rf_samples = [float(raw(x)[4]["rf_mhz"][0]) for x in angles]
                trial_angles = list(angles)
                # Include exact RF-boundary designs; useful for both fixed and ranged RF.
                for target in set([s.rf_min_mhz, s.rf_max_mhz]):
                    for j in range(len(angles)-1):
                        if (rf_samples[j]-target)*(rf_samples[j+1]-target) > 0:
                            continue
                        try:
                            root_angle = brentq(lambda x: float(raw(x)[4]["rf_mhz"][0])-target,
                                                angles[j], angles[j+1], xtol=1e-10)
                        except ValueError:
                            continue
                        trial_angles.append(root_angle)
                family_candidates = []
                for x in trial_angles:
                    rec = record(x)
                    if rec:
                        candidates.append(rec)
                        family_candidates.append((rec["gain_per_rf_w"], x))
                        scan.append(dict(tilt_deg=rec["acoustic_tilt_from_basal_deg"], mode=mode,
                                         interaction=pi+rec["diffracted_polarization"],
                                         gain_per_rf_w=rec["gain_per_rf_w"], m2_s3_kg=rec["m2_s3_kg"], rf_mhz=rec["rf_mhz"]))
                if not fixed:
                    # Constrained feasible objective and several local intervals.
                    def objective(x):
                        rec = record(x)
                        return -rec["gain_per_rf_w"] if rec else 1e9
                    for _, x in sorted(family_candidates, reverse=True)[:8]:
                        left, right = max(lo, x-2*s.grid_deg), min(hi, x+2*s.grid_deg)
                        result = minimize_scalar(objective, bounds=(left, right), method="bounded",
                                                 options={"xatol": 1e-9})
                        if result.success:
                            rec = record(result.x)
                            if rec:
                                candidates.append(rec)
    return candidates, scan


def build_report(settings, model, candidates, scan):
    if not candidates:
        raise ValueError("No admissible design found. Broaden RF/cut bounds or reduce grid spacing; "
                         "the parallel-tangent constraint can also make fixed RF infeasible.")
    candidates.sort(key=lambda r: r["gain_per_rf_w"], reverse=True)
    unique = []
    for r in candidates:
        if all((r["acoustic_mode_index"], r["incident_polarization"], r["diffracted_polarization"]) !=
               (u["acoustic_mode_index"], u["incident_polarization"], u["diffracted_polarization"]) or
               abs(r["acoustic_tilt_from_basal_deg"]-u["acoustic_tilt_from_basal_deg"]) > 0.5 or
               abs(r["incident_wavevector_signed_to_c_deg"]-u["incident_wavevector_signed_to_c_deg"]) > 0.5
               for u in unique):
            unique.append(r)
        if len(unique) == 20:
            break
    warnings = [
        "Chang 1976 review text audited; modern numerical extensions are documented. Historical device data are not reproduced.",
        "One specified interaction plane only; no certified global optimum across all crystal orientations.",
        "Optical activity omitted; both optical wavevectors restricted to at least 5 degrees from [001].",
        "Photoelastic coefficients held at their 632.8-nm values; wavelength-dependent efficiency is approximate.",
        "Uniform two-wave Bragg coupling assumed; finite beams, acoustic diffraction, and full walkoff overlap omitted.",
        "Face refraction, transducer impedance/resonance, heating and coatings are not designed by this program.",
    ]
    if not settings.strain_only:
        warnings.append("Rotation term assumes photoelastic inputs are intrinsic strain coefficients; effective measured coefficients may already include rotation.")
    if unique[0]["klein_cook_q_diagnostic"] < 10:
        warnings.append("Scalar Klein-Cook Q diagnostic is below 10; the two-wave approximation needs a multi-order check.")
    if settings.loss_db_per_mm_100mhz == 0:
        warnings.append("Acoustic attenuation is assumed zero; supply measured loss and path length for a practical estimate.")
    if settings.optical_transmission == 1:
        warnings.append("Optical transmission is assumed unity; throughput excludes unprovided surface and absorption losses.")
    if settings.wavelength_nm > 1000:
        warnings.append("Sellmeier dispersion extrapolated beyond the verified 400-1000-nm fit interval; absorption not modeled.")
    return dict(version=VERSION, scope=SCOPE, inputs=asdict(settings), material=asdict(model.m),
                literature_audit="See AUDIT.md: Chang 1976 Sections II and V; explicit rotation and flux extensions.",
                warnings=warnings, optimum=unique[0], alternatives=unique[1:], coarse_cut_scan=scan)


def write_outputs(report, prefix):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(str(prefix)+".json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    rows = [report["optimum"]]+report["alternatives"]
    fields = [k for k, v in rows[0].items() if not isinstance(v, list)]
    with Path(str(prefix)+".csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    scan = report["coarse_cut_scan"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
    keys = sorted(set((r["mode"], r["interaction"]) for r in scan))
    for mode, interaction in keys:
        ss = [r for r in scan if (r["mode"], r["interaction"]) == (mode, interaction)]
        best = {}
        for r in ss:
            if r["tilt_deg"] not in best or r["gain_per_rf_w"] > best[r["tilt_deg"]]["gain_per_rf_w"]:
                best[r["tilt_deg"]] = r
        ss = [best[x] for x in sorted(best)]
        x = [r["tilt_deg"] for r in ss]
        axs[0].plot(x, [r["gain_per_rf_w"] for r in ss], label=f"mode {mode}, {interaction[0]}→{interaction[1]}")
        axs[1].plot(x, [r["rf_mhz"] for r in ss])
    for ax in axs:
        ax.axvline(report["optimum"]["acoustic_tilt_from_basal_deg"], color="black", ls="--", lw=1)
        ax.set_xlabel("Acoustic normal tilt from basal plane (degrees)")
        ax.grid(alpha=0.25)
    axs[0].set_ylabel("First-lobe coupling per RF watt (1/W)")
    axs[0].legend(fontsize=8)
    axs[1].set_ylabel("Phase-matched RF frequency (MHz)")
    fig.suptitle(f"TeO₂ {report['inputs']['wavelength_nm']:g} nm — {report['inputs']['design_mode']} cut search (v{VERSION})")
    fig.tight_layout()
    fig.savefig(str(prefix)+".png", dpi=170)
    plt.close(fig)
    write_response_graphs(report, prefix)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--wavelength-nm", type=float, required=True)
    p.add_argument("--rf-mhz", type=float, help="Fix RF frequency; overrides RF bounds")
    for name, default in asdict(Settings(632.8)).items():
        if name == "wavelength_nm":
            continue
        if isinstance(default, bool):
            p.add_argument("--"+name.replace("_", "-"), action="store_true", default=default)
        elif isinstance(default, str):
            choices = {"interaction": ["anisotropic", "isotropic", "all"],
                       "incident": ["o", "e", "both"], "acoustic_mode": ["slow", "middle", "fast", "all"],
                       "design_mode": ["efficiency", "wide-aperture"]}
            p.add_argument("--"+name.replace("_", "-"), choices=choices[name], default=default)
        else:
            p.add_argument("--"+name.replace("_", "-"), type=float, default=default)
    p.add_argument("--material-json", type=Path, help="Override room-temperature scalar material constants")
    p.add_argument("--output-prefix", default="teo2_result")
    args = p.parse_args()
    raw = vars(args).copy()
    prefix, material_path, rf = (raw.pop(k) for k in ("output_prefix", "material_json", "rf_mhz"))
    if rf is not None:
        raw["rf_min_mhz"] = raw["rf_max_mhz"] = rf
    try:
        material = Material(**json.loads(material_path.read_text())) if material_path else Material()
        report = optimize(Settings(**raw), material)
        write_outputs(report, prefix)
    except (ValueError, TypeError, OSError) as e:
        p.error(str(e))
    r = report["optimum"]
    print(f"TeO2 {args.wavelength_nm:g} nm: {SCOPE}")
    print(f"Acoustic normal to C-axis: {r['acoustic_normal_to_c_axis_deg']:.3f} deg")
    print(f"Signed tilt from basal plane: {r['acoustic_tilt_from_basal_deg']:.3f} deg")
    print(f"Optical k to C-axis: {r['incident_wavevector_to_c_axis_deg']:.3f} deg")
    print(f"Polarization: {r['incident_polarization']} -> {r['diffracted_polarization']}; acoustic mode {r['acoustic_mode_index']}")
    print(f"RF: {r['rf_mhz']:.4f} MHz; phase speed: {r['acoustic_phase_speed_m_s']:.2f} m/s")
    print(f"M2: {r['m2_s3_kg']/1e-15:.2f} x 10^-15 s^3/kg")
    print(f"First peak RF power: {r['rf_first_peak_w']:.5f} W")
    print(f"Selected RF power: {r['rf_selected_w']:.5f} W (within budget)")
    print(f"Estimated intrinsic efficiency: {100*r['intrinsic_efficiency']:.2f}%")
    print(f"Estimated throughput efficiency: {100*r['throughput_efficiency']:.2f}%")
    print(f"Acoustic walkoff: {r['acoustic_walkoff_deg']:.3f} deg")
    print(f"Saved {prefix}.json, {prefix}.csv, cut graph, four performance graphs and sweep CSVs")
    for warning in report["warnings"]:
        print("MODEL LIMIT: "+warning)


if __name__ == "__main__":
    main()
