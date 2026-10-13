# Validation — version 0.2.0, 12 October 2026

Fourteen numerical checks passed. This verifies mathematical identities and numerical behavior, not experimental device performance or reproduction of a historical Chang device.

## Executed examples

All examples assume 5 mm boundary-normal interaction thickness, 5 mm × 1 mm acoustic area, 0.5 conversion, unit intensity overlap, 0.05 W RF budget, 0.95 optical transmission and zero acoustic attenuation.

| Quantity | 532 nm / 80 MHz, efficiency | 532 nm / 80 MHz, wide aperture | 632.8 nm / ranged RF |
|---|---:|---:|---:|
| Acute face-plane angle to C-axis (degrees) | 4.1895 | 5.7456 | 1.3939 |
| Incident internal wavevector angle (degrees) | 20.3927 | 13.2490 | 10.7785 |
| RF carrier (MHz) | 80.0000 | 80.0000 | 30.5551 |
| Acoustic phase speed (m/s) | 635.35 | 648.81 | 621.41 |
| M2 (10^-15 s^3/kg) | 775.01 | 760.48 | 738.29 |
| First-peak RF input (W) | 0.07301 | 0.07440 | 0.10847 |
| Intrinsic efficiency (%) | 92.84 | 92.20 | 76.65 |
| Throughput (%) | 88.20 | 87.59 | 72.82 |
| Acoustic walkoff (degrees) | 33.661 | 41.122 | 13.080 |
| Parallel-tangent residual | -1.378e-02 | -2.159e-16 | 9.054e-03 |
| Optical flux correction | 1.000391 | 1.000400 | 1.000092 |

Reciprocal and mirrored orientations can have the same coupling. The directed acoustic normal, signed incident angle and polarization are preserved in each JSON report; the acute face-plane angle alone is insufficient for manufacturing. Large acoustic walkoff must be checked against actual finite-beam overlap.

## Convergence

The 532 nm efficiency example was repeated with 1-degree and 2-degree initial grids. Comparing mirrored cuts by absolute tilt magnitude:

| Quantity | Absolute difference |
|---|---:|
| Acoustic tilt magnitude (degrees) | 5.847e-07 |
| M2 (s^3/kg) | 1.975e-24 |
| First-peak RF power (W) | 1.832e-15 |
| Throughput efficiency fraction | 7.883e-15 |

This demonstrates convergence for the example, not physical accuracy at the printed precision or global optimality outside the selected design domain.

## Checks completed

- Positive elastic stiffness and acoustic eigenproblem residuals.
- Slow-shear speed scale and engineering shear convention.
- Optical surface limits and reproducible fit to primary measured indices.
- Exact geometric momentum matching and figure-of-merit scale.
- Nominal detuned response agrees with matched efficiency.
- Zero power, first peak, overdrive zero and low-power slope.
- RF detuning and near-axis masking.
- Finite rotated-tensor derivative agrees with the explicit rotation term.
- Rotation vanishes in an isotropic optical tensor.
- Poynting flux reduces to cos(theta)/n^3 for isotropic oblique waves.
- Closed-form efficiency agrees with independent matrix-exponential evolution of a lossless two-mode system.
- Parallel-tangent geometry has zero first-order angular mismatch, with nonzero second-order mismatch.
- Fixed-frequency wide-aperture optimizer satisfies both RF and tangent constraints.
- Invalid numerical inputs are rejected.

The fourteen tests include multiple related checks within each test method.

## Source and physical limits

The 1976 review text was audited, but unclear scanned angle factors were handled through an explicit independent flux derivation. See AUDIT.md. Optical activity, photoelastic dispersion and finite-beam field overlap remain omitted. The rotational contribution depends on whether supplied coefficients are intrinsic or already effective. No experimental efficiency error bars are available.
