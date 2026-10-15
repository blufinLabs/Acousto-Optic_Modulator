# Validation — version 0.2.0, 12 October 2026

## RF input-match requirement update — 13 October 2026

The current requirement is inclusive: input return loss >=30 dB (S11 <=-30 dB),
reflected power <=0.1%, and VSWR <=1.065:1 over [f_min, f_max], with port 2
terminated in 50 ohms. The exact VSWR corresponding to 30 dB is 1.0653109:1;
the specified 1.065:1 limit is slightly stricter. Numerical band endpoints
remain to be agreed. A supplied scalar VSWR does not establish band compliance.
This is an acceptance requirement, not a measured result or optical insertion
loss.

## Engineering-package extension — 13 October 2026

Validation here concerns numerical behavior and generated artifacts, **not
manufacturing validation**. The original fourteen physics checks passed after
installing the existing `requirements.txt` dependencies. Fifteen physics checks
then passed with the new winning-azimuth/old-report compatibility test.
Package-search, geometry, CLI and artifact test results are recorded below
after execution.

The first integrated 34-test run reported 1 failure and 3 errors: the initial
5 mm package interaction fixture was not contained by the bounded acoustic
field in the assumed 10 mm crystal. A 3 mm slab with 7 × 3 mm acoustic area
was verified by the package-layer checks and adopted as the package-only
default. The legacy optimizer's 5 mm default and explicit user dimensions
remain unchanged.

Representative package-generation smoke command (example, not actual user
operating values; outputs intentionally outside the repository):

```sh
python teo2_optimizer.py --design-package --wavelength-nm 532 --rf-mhz 80 \
  --design-mode wide-aperture --acoustic-mode slow --incident e --grid-deg 5 \
  --azimuth-min-deg 30 --azimuth-max-deg 60 --azimuth-step-deg 15 \
  --rf-ceiling-w 0.9 --duty-cycle-percent 100 --output-prefix /tmp/teo2_example
python -m unittest discover -v
```

RF perturbations are ±0.1% at the selected power and **fixed** cut/alignment.
Nominal RF remains a hard equality for the winner and alternatives. Tolerance
response evaluations are explicitly conditional bulk predictions, not measured
yield or a coupled worst-case tolerance analysis. Finite Gaussian-envelope
containment is not full overlap validation; infinite tails are excluded according
to the reported threshold.

**Historical requirement:** passive VSWR is Vmax/Vmin ≥1. Reversal from 1000:1 to
1:1000 is not better or valid; the previous suggestion was wrong. The user's
subsequent “insertion loss (VSWR) >30 dB” request supersedes the prior
provisional ≤1.5:1 target. The user clarified **less than 0.1% of forward RF
power reflected**, confirming **RF return loss >30 dB**, equivalent to
**VSWR <1.0653109:1**. The requirement is strict, not ≥30 dB. It is not optical
insertion loss or a measured result. Nominal 50 Ω alone cannot establish
reflected or net accepted RF power. The 1 W illustration does not relax the
strictly sub-1-W power ceiling.
The latest requirement above supersedes this historical strict-inequality
interpretation.

No supplier capability, cold-weld process, LiNbO₃ conversion, coating,
temperature rise, continuous-duty endurance, switching time or extinction was
qualified. Supplier-reference access failed DNS resolution. Numerical fallback
tolerances remain engineering proposals pending independent primary-source
verification and Crys-Teh review. Generated packages remain **NOT FOR FABRICATION**.

## Original v0.2 validation

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
