# Chang review audit and v0.2 model changes

## What was obtained

Readable text of the original 1976 IEEE review was obtained from a scanned copy.
Sections II and V were examined. The 1974 paper's abstract and the 1975 publication
record were found; complete texts of those papers and the 1977 overview were not
obtained. Therefore this package no longer states that all Chang text was unavailable,
but it also does not claim complete historical-paper or historical-device validation.

## Audited model relationships

| Item | Finding and action |
|---|---|
| Two-wave response | The implemented detuned sine-squared response has the functional structure of the review's Eq.13. Its numerical implementation is independently checked against evolution of a lossless 2x2 coupled-amplitude system. |
| Momentum mismatch | The review describes mismatch normal to the interaction boundary. v0.1 rotated that coordinate with the input alignment sweep. v0.2 fixes the boundary at the nominal design direction throughout each sweep. |
| Angular factors | Some scanned coupling factors are unclear in OCR. v0.2 derives explicit flux factors from E cross H rather than claiming the OCR was transcribed accurately. |
| Rotation | The review cautions about rotational contributions for shear waves in birefringent media. v0.2 adds a separately reported kinematic rotation term under an explicit intrinsic-coefficient assumption. |
| Wide aperture | The parallel-tangent construction is imposed in the optional wide-aperture mode. Efficiency optimization remains a separate objective. |
| Higher orders | The two-wave approximation remains an assumption. A scalar Klein-Cook diagnostic is added, with a low-Q warning. |

## Derivation of the added rotation term

For displacement gradient H=u outer q, separate symmetric strain S=(H+H.T)/2
from antisymmetric rotation W=(H-H.T)/2. An infinitesimal physical rotation has
R=I+epsilon*W, so the laboratory impermeability tensor transforms as

`R*B0*R.T = B0 + epsilon*(W*B0-B0*W) + O(epsilon^2)`.

This term is symmetric and vanishes when B0 is isotropic. A centered finite
rotation with a matrix exponential reproduces the implementation in a regression
check. The derivation is a modern mathematical extension; it is not a numerical
formula copied from a Chang table. Existing effective coefficients may already
contain a rotational contribution. The program offers `--strain-only` and reports
both contributions rather than concealing that uncertainty.

## Derivation of optical boundary flux

For a unit displacement eigenvector D, set E=B0*D. Maxwell's plane-wave relation
gives H proportional to n*s cross E. The energy flux projection through a boundary
normal z is proportional to

`J=n*((s.z)*|E|^2-(s.E)*(E.z))`.

Normalization of the two mode amplitudes by their boundary flux gives a symmetric
coupling factor proportional to `p_eff/sqrt(J_i*J_d)`. Consequently

`kappa^2=(pi/lambda)^2*I_a*p_eff^2/(2*rho*v^3*J_i*J_d)`.

For an isotropic index, J=cos(theta)/n^3. At normal incidence this recovers the
usual `n^6*p_eff^2/(rho*v^3)` material factor. This normalization is checked in
the isotropic oblique limit. The selected interaction thickness is normal to the
reference boundary; it is not the distance along a moving energy ray.

## Exact wide-aperture construction

Normals to the ordinary sphere and extraordinary ellipsoid are parallel when

`tan(theta_e)=(n_e/n_o)^2*tan(theta_o)`.

Their wavevector difference determines the acoustic normal. The program then
calculates the correct phase velocity and phase-matched RF carrier. A finite
alignment derivative check verifies that first-order longitudinal mismatch vanishes
at a returned parallel-tangent design. The remaining mismatch is second order;
finite aperture and finite beam overlap still matter.

## What remains unverified

Absolute efficiency has not been tested against a fabricated device or Chang's
historical examples. Optical activity, effective-coefficient calibration, acoustic
field overlap and multi-order coupling remain material limits. Historical approximate
cut formulas are not silently substituted for the exact numerical construction.

## Engineering-package extension

The original material data and bulk coupling equations are retained. Package
mode adds a separate engineering-settings layer, fixed-RF azimuth-family search,
closed-solid geometry and reporting, not a validated transducer or coating model.
Azimuth varies the common coplanar optical/acoustic family, not three independent
wavevectors. Candidate filtering precedes report deduplication; refined cuts and
alternatives are rechecked for nominal RF, original optical constraints and
physical package feasibility. The selected azimuth is carried into response
sweeps, while old reports without candidate-level azimuth use input azimuth.
Cut-scan curves are grouped by azimuth instead of conflating different planes.

### Boundary direction calculation

At a translationally invariant planar dielectric boundary, phase continuity
conserves tangential wavevector. For an internal dimensionless vector
`b = k λ/(2π)` and outward exit normal `n`, define
`bt = b - (b·n)n`. In air (index approximated as 1), a propagating output requires
`|bt|² < 1`; the outgoing vector is
`b_air = bt + sqrt(1-|bt|²)n`. The internal extraordinary dispersion surface and
its gradient supply the internal wavevector and energy ray separately. This
construction does not assume the internal extraordinary energy ray satisfies
scalar Snell's law. TIR/evanescent cases and grazing output are non-transmitting.
It is a kinematic boundary calculation; Fresnel amplitudes, coating behavior,
air-index corrections and optical activity remain outside the model.

### Specification provenance and release audit

Official Crys-Teh, Newport and Edmund reference retrieval failed DNS resolution
during this task. Consequently no new verified supplier numerical specifications
are claimed. Fallback tolerance numbers are proposed engineering choices with
rationale, **not sourced vendor capabilities or universal standard tolerances**.
Primary-source verification and material-specific supplier acceptance are
explicit release blockers. Standards notation does not establish attainable
TeO₂ cut, polish or coating quality. Coating validity must be established at
the user wavelength, both face incidence angles and selected polarization.

Requested LiNbO₃/vacuum cold-weld technology does not substantiate a cut,
thickness, metallurgy, process, impedance or conversion efficiency. Legacy
conversion/overlap/transmission remain conditional model inputs. Average power
is not a thermal calculation. Generated concept drawings are NOT FOR FABRICATION;
no generation, convergence or boolean approval can qualify the hardware.

Passive VSWR is Vmax/Vmin ≥1. Reversing 1000:1 to 1:1000 is invalid, not improved
matching. User intent was unresolved; ≤1.5:1 was only a provisional proposed target.
Nominal 50 Ω is not evidence of measured VSWR. Linewidth, switching targets,
extinction and temperature requirements remain unspecified.

The subsequent user request, “insertion loss (VSWR) >30 dB,” superseded the
old provisional target. VSWR, RF return loss and optical insertion loss are
different metrics. The user then confirmed “for 1 watt inserted, less than 0.1%
is returned,” which was recorded as strict >30 dB return loss. The current
requirement supersedes that strict inequality: input return loss >=30 dB
(S11 <=-30 dB), reflected power <=0.1%, and VSWR <=1.065:1 across [f_min,
f_max], with port 2 terminated in 50 ohms. Return loss is
`-20 log10(|Gamma|)`; exactly 30 dB corresponds to VSWR 1.0653109:1, so the
stated 1.065:1 cap is slightly stricter. The numerical band endpoints remain
to be agreed. None of these requirements is optical insertion loss or measured
RF-match evidence. The 1 W example does not change the strictly sub-1-W
operating power ceiling.
