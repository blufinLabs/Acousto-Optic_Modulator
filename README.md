# TeO₂ acousto-optic optimizer — version 0.2.0

This Python program selects a crystal orientation for a supplied optical wavelength
and predicts first-order diffraction efficiency using a uniform two-wave Bragg model.
It produces the four requested graphs: efficiency versus wavelength, RF carrier
frequency, internal angular alignment, and RF input power.

The updated model uses the coupled-wave framework described in I. C. Chang's
1976 review, adds exact numerical optical phase matching and acoustic eigenmodes,
and offers a parallel-tangent wide-aperture option. The source audit distinguishes
what was read in Chang from the mathematical extensions implemented here.
**This is an engineering prediction, not experimentally validated device performance.**

## Install and run

Python 3.10 or newer is required. Extract the package, open PowerShell in the
`teo2_optimizer` folder, and run:

```powershell
py -m pip install -r requirements.txt
py teo2_optimizer.py --wavelength-nm 532 --rf-mhz 80 --rf-budget-w 0.05 --optical-transmission 0.95
```

On macOS/Linux substitute `python3` for `py`.

Maximum coupling per RF watt is the default objective. To impose Chang's
parallel-tangent geometry for wide angular acceptance:

```powershell
py teo2_optimizer.py --wavelength-nm 532 --rf-mhz 80 --design-mode wide-aperture --rf-budget-w 0.05 --optical-transmission 0.95 --output-prefix wide_design
```

To let the optimizer choose the phase-matched RF frequency:

```powershell
py teo2_optimizer.py --wavelength-nm 632.8 --rf-min-mhz 20 --rf-max-mhz 200 --output-prefix red_design
```

For a comparison using only the strain contribution to the photoelastic coupling:

```powershell
py teo2_optimizer.py --wavelength-nm 532 --rf-mhz 80 --strain-only --output-prefix strain_comparison
```

The rotation contribution is included by default. Use `--strain-only` if supplied
calibrated effective coefficients already incorporate rotation, or to quantify
sensitivity to that modeling assumption. v0.2 strain-only results still use the
new boundary and optical-flux treatment; they are not identical to v0.1 results.

```powershell
py teo2_optimizer.py --help
py -m unittest discover -v
```

No internet access is needed after installing dependencies.

## Inputs and assumptions

| Input | Default / meaning |
|---|---|
| `--wavelength-nm` | Required; vacuum optical wavelength |
| `--rf-mhz` | Optional fixed RF carrier; overrides the RF range |
| `--rf-min-mhz`, `--rf-max-mhz` | 20–200 MHz when fixed RF is absent |
| `--design-mode` | `efficiency` or `wide-aperture` |
| `--length-mm` | 5 mm interaction slab thickness along its boundary normal |
| `--acoustic-width-mm`, `--acoustic-height-mm` | 5 mm × 1 mm acoustic wavefront/transducer area |
| `--rf-budget-w` | 1 W available input; RF power can be reduced to the first efficiency peak |
| `--conversion` | 0.5 RF-to-acoustic power conversion, assumed constant over a sweep |
| `--overlap` | 1.0 acoustic intensity factor; not a calculated spatial-overlap integral |
| `--optical-transmission` | 1.0; includes surface/absorption loss only when supplied |
| `--loss-db-per-mm-100mhz` | 0; user-supplied intensity loss at 100 MHz, scaled as RF frequency squared |
| `--acoustic-path-mm` | 0; path used with the attenuation coefficient |
| `--max-walkoff-deg` | 85° acoustic energy walkoff limit |
| `--azimuth-deg` | 45°, the [110]–[001] interaction plane |
| `--tilt-min-deg`, `--tilt-max-deg` | -30° to +30° acoustic wave-normal tilt from the basal direction |
| `--theta-min-deg`, `--theta-max-deg` | 10° to 85° internal optical wavevector angles to the C-axis |
| `--incident` | `both`, `o`, or `e` |
| `--interaction` | `anisotropic`, `isotropic`, or `all`; wide-aperture requires anisotropic |
| `--acoustic-mode` | `all`, `slow`, `middle`, or `fast`, ordered by phase speed |
| `--grid-deg` | 2° coarse starting grid; repeat with 1° to assess convergence |
| `--material-json` | Optional overrides of the scalar `Material` constants |

The example conversion, dimensions, power budget and transmission are assumptions,
not measurements. RF carrier frequency is distinct from RF-envelope modulation bandwidth.
The default unity transmission and overlap can produce ideal values near 100%.

Built-in optical dispersion is fitted over 400–1000 nm. Exploratory extrapolation
up to 4500 nm requires `--allow-dispersion-extrapolation`; infrared absorption and
photoelastic dispersion are not provided. For example:

```powershell
py teo2_optimizer.py --wavelength-nm 1064 --allow-dispersion-extrapolation --rf-mhz 80
```

## Crystal orientation and angle definitions

Crystal Cartesian axes are X=[100], Y=[010], Z=[001], with Z the C-axis.
For basal azimuth phi, define `t=(cos(phi),sin(phi),0)`.
The acoustic normal is `q=cos(alpha)*t+sin(alpha)*Z`, and the optical wavevector
is `s=sin(theta)*t+cos(theta)*Z`. Both signs of theta are searched.

| Report field | Interpretation |
|---|---|
| `acoustic_normal_to_c_axis_deg` | Directed angle from q to +[001]: 90-alpha |
| `acoustic_face_plane_to_c_axis_deg` | Acute angle between transducer face plane and C-axis: abs(alpha) |
| `acoustic_tilt_from_basal_deg` | Signed acoustic normal tilt from t toward +[001] |
| `incident_wavevector_to_c_axis_deg` | Internal incident wavevector angle to +[001] |
| `diffracted_wavevector_to_c_axis_deg` | Internal diffracted wavevector angle to +[001] |
| `acoustic_normal_xyz` | Unit face normal, removing angle-convention ambiguity |
| `optical_ray_separation_deg` | Separation of optical energy rays, which differ from wavevector directions |

The transducer face is perpendicular to q. The interaction boundary normal is
chosen along the nominal incident wavevector for each candidate, then kept fixed
through its response sweeps. Optical entrance/exit refraction is not calculated.
These are internal angles, not external alignment angles or a fabrication drawing.

The optimization is within one selected crystal plane and the supplied bounds.
It is not a certified global search of three-dimensional crystal orientation.
Numerical angle precision is not manufacturing accuracy. Reciprocal polarization
geometries and mirrored directions can have equally good objectives.

## Two design modes

`efficiency` maximizes coupling per RF watt subject to the selected RF, mode and
orientation bounds. With RF power adjustable from zero to the budget, this also
maximizes attainable efficiency. If several cuts can reach the ideal peak,
it chooses the one requiring the least power. The report separately gives
`efficiency_at_exact_budget`, exposing overdrive at the full budget.

`wide-aperture` imposes parallel tangents to the ordinary and extraordinary
wavevector surfaces, then maximizes coupling among that family of designs.
This makes first-order angular phase mismatch vanish at the design point;
it does not guarantee unlimited aperture or a specified finite-beam efficiency.
In the selected plane, the exact relation is

`tan(theta_e) = (n_e/n_o)^2 * tan(theta_o)`.

The acoustic vector comes from the difference between the two optical wavevectors.
The program solves the required RF frequency using the actual acoustic eigenmode
speed. When RF is fixed, it solves the one-parameter family for matching roots.
When RF is ranged, it compares feasible points and refines candidate intervals.
It uses full non-gyrotropic index surfaces, rather than a small-birefringence cut-angle
approximation. If no feasible parallel-tangent design exists, it reports an error
instead of silently falling back to the efficiency mode.

## Equations implemented

All internal quantities use SI units. Optical frequency shift is neglected when
computing indices because RF is tiny compared with optical frequency.

1. **Acoustics:** solve the Christoffel eigenproblem for phase speed v and unit
   displacement u. Group velocity comes from the phase-speed gradient and gives
   the acoustic walkoff diagnostic.
2. **Optical indices:** `indices()` uses an independently reproducible weighted
   two-oscillator fit to Uchida's measured Table I. `fit_dispersion.py` includes
   those measured data and reproduces the fit. The 9.24 and 4.70 eV resonances
   are fixed. This is our refit, not a transcription of Uchida's original parameter table.
3. **Exact geometric phase matching:** with `b=k*lambda/(2*pi)`, solve
   `b_d=b_i+Q*q` on the diffracted wavevector surface. The ordinary surface is
   a sphere; the extraordinary surface satisfies
   `b_x^2/n_e^2+b_y^2/n_e^2+b_z^2/n_o^2=1`.
   RF is `abs(Q)*v/lambda`. Both frequency-shift signs are compared.
4. **Strain coupling:** `S=(u outer q+q outer u)/2`. The photoelastic matrix maps
   `[Sxx,Syy,Szz,2Syz,2Sxz,2Sxy]` into `[Bxx,Byy,Bzz,Byz,Bxz,Bxy]`.
5. **Crystal rotation:** define `W=(u outer q-q outer u)/2`,
   `B0=diag(1/n_o^2,1/n_o^2,1/n_e^2)`, and
   `delta_B_rotation=W*B0-B0*W`. This follows by linearizing `R*B0*R.T`.
   It vanishes for an optically isotropic material. The effective coefficient is
   `p_eff=D_d.T*(delta_B_strain+delta_B_rotation)*D_i` for unit displacement
   eigenvectors D. The strain and rotation terms are reported separately.
6. **Optical flux:** for each unit D, `E=B0*D`, and the normalized flux through
   the fixed boundary normal z is
   `J=n*((s.z)*|E|^2-(s.E)*(E.z))`.
   In an isotropic medium this reduces to `cos(theta)/n^3`.
7. **Coupling:** acoustic intensity is
   `I_a=P_RF*conversion*overlap*10^(-loss_dB/10)/area`.
   With boundary-normal interaction thickness L,
   `kappa^2=(pi/lambda)^2*I_a*p_eff^2/(2*rho*v^3*J_i*J_d)`.
   The conventional diagnostic `M2=n_i^3*n_d^3*p_eff^2/(rho*v^3)` is reported
   separately from its geometric flux correction.
8. **Matched efficiency:** `eta_internal=sin(kappa*L)^2` and
   `eta_throughput=T*eta_internal`. Power is reduced to the first peak if the
   budget is larger. T is the supplied optical transmission.
9. **Detuned response:** conserve transverse momentum at the fixed boundary,
   solve the diffracted surface for its forward longitudinal component, and
   calculate longitudinal mismatch Delta k. Then
   `eta=T*kappa^2/(kappa^2+(Delta k/2)^2)`
   `       *sin(L*sqrt(kappa^2+(Delta k/2)^2))^2`.

The two-wave response has the functional form of Chang's 1976 Eq.13. Rotation
and explicit Poynting normalization are documented mathematical extensions;
we do not claim they were numerically transcribed from the review's scanned Eq.14.
See `AUDIT.md` for the derivation checks and remaining uncertainty.

## Outputs and graphs

Each run writes:

- `<prefix>.json`: inputs, material constants, optimum, alternatives, model limits.
- `<prefix>.csv`: numerical candidate designs.
- `<prefix>.png`: coupling per RF watt and phase-matched frequency versus cut tilt.
- `<prefix>_performance.png` and `.pdf`: four-panel performance figure.
- Four sweep CSV files with independent variable, efficiency fraction and mismatch.

| Curve | Parameters held fixed |
|---|---|
| Efficiency versus wavelength | Cut, RF, internal direction, RF power |
| Efficiency versus RF frequency | Cut, wavelength, internal direction, RF power |
| Efficiency versus angular alignment | Cut, wavelength, RF, RF power, physical interaction boundary |
| Efficiency versus RF power | Cut, wavelength, RF and internal direction |

These are fixed-design sensitivity curves. Each sample is not reoptimized.
Angular alignment means an internal wavevector rotation in the selected plane.
The power curve extends beyond the first peak to illustrate overdrive; it is not
an RF power rating. Wavelength spans +/-1.5%, RF +/-10%, and angle +/-1 degree.
Code-level `response()` supports custom sweeps. Unmodeled near-axis samples
appear as NaN gaps, not predictions of zero efficiency.

## Material data and references

The default room-temperature data set is from Mys et al. (2019), section 2:
density 5990 kg/m^3; `(C11,C12,C13,C33,C44,C66)` =
`(53.2,48.6,21.2,108.5,24.4,55.2)` GPa; and
`(p11,p12,p13,p31,p33,p44,p66)` =
`(0.0074,0.187,0.340,0.0905,0.240,-0.170,-0.0463)`.
Photoelastic coefficients are held at their 632.8 nm values.

| Reference | Access / use |
|---|---|
| Chang, *Noncollinear acousto-optic filter with large angular aperture*, APL 25, 370–372 (1974), [DOI](https://doi.org/10.1063/1.1655512) | Abstract accessed; wide-aperture historical reference, not a full numerical reproduction. |
| Chang, *Analysis of the noncollinear acousto-optic filter*, Electronics Letters 11, 617–618 (1975), [DOI](https://doi.org/10.1049/el:19750470) | Publication record/abstract; complete equations not obtained. |
| Chang, *Acoustooptic Devices and Applications*, IEEE Trans. Sonics Ultrasonics SU-23(1) (1976), [publisher record](https://ieeexplore.ieee.org/document/1538911), [readable scanned-paper text](https://www.scribd.com/document/71916308/Acoustooptic-Devices-and-Applications) | Primary article text examined: Section II two-wave efficiency and phase mismatch; Section V noncollinear wide-aperture geometry. Some equation typography is unreliable in OCR. |
| Chang, *Tunable Acousto-Optic Filters: An Overview*, Optical Engineering 16, 455–460 (1977), [DOI](https://doi.org/10.1117/12.7972116) | Bibliographic reference; full text unavailable. |
| Uchida, Phys. Rev. B 4, 3736–3745 (1971), [DOI](https://doi.org/10.1103/PhysRevB.4.3736), [author-uploaded paper](https://www.researchgate.net/publication/235448791_Optical_Properties_of_Single-Crystal_Paratellurite_TeO2) | Primary measured index table; refitted independently. |
| Mys et al., Ukr. J. Phys. Opt. 20 (2019), [publisher PDF](https://ifo.lviv.ua/journal/UJPO_PDF/2019_1/0401_2019.pdf) | Material constants and acoustic eigenmode formulation. |
| Mys et al., Ukr. J. Phys. Opt. 15 (2014), [publisher PDF](https://ifo.lviv.ua/journal/UJPO_PDF/2014_3/0403_2014.pdf) | Literature velocity and figure-of-merit scale checks. |
| Gass and Sambles, *Accurate design of a noncollinear acousto-optic tunable filter*, Optics Letters 16, 429–431 (1991), [DOI](https://doi.org/10.1364/OL.16.000429), [paper text](https://www.researchgate.net/publication/26831209_Accurate_design_of_a_noncollinear_acousto-optic_tunable_filter) | Exact parallel-tangent geometry cross-check; documents limitations of approximate birefringence designs. |

No source PDFs are redistributed.

## Practical limits and validation

Optical activity/gyrotropy is still omitted. Both optical beams are excluded within
5 degrees of [001], with a 10-degree default search bound. This is a modeling rule,
not a universal measured accuracy boundary. Near-axis elliptical/circular eigenstates
require a gyrotropic model and may change the optimum.

Measured effective photoelastic coefficients can mix intrinsic strain and rotation;
adding an explicit rotation term to such values can double-count it. The default
interpretation and `--strain-only` comparison are therefore explicitly reported.
No experimental data establish a quantitative efficiency error bar here.

Finite beam averaging, acoustic diffraction, full spatial overlap, transducer impedance
and resonance, heating, and external face refraction are not modeled. Large acoustic
walkoff can make the assumed interaction length physically unrealizable.
The scalar Klein-Cook Q is a diagnostic for the two-wave approximation, not a
sufficient multi-order validity test for anisotropic diffraction. Low values emit a warning.
The characteristic transit-time diagnostic is not a measured optical rise-time prediction.

Fourteen numerical checks pass, including independent coupled-amplitude evolution,
finite tensor rotation, Poynting-flux limits, exact phase matching, primary-index fit,
and angular stationarity of the parallel-tangent design. The examples include both
532 nm/80 MHz design modes and a 632.8 nm ranged-RF design. `VALIDATION.md` contains
executed results and convergence checks; `AUDIT.md` explains the changes from v0.1.
