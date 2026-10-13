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
