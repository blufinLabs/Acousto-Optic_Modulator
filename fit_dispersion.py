"""Reproduce the program's optical-index fit from PRIMARY measured data.

Uchida, Phys. Rev. B 4, 3736 (1971), Table I: observed indices at 20 C.
https://doi.org/10.1103/PhysRevB.4.3736
Author-uploaded full paper:
https://www.researchgate.net/publication/235448791_Optical_Properties_of_Single-Crystal_Paratellurite_TeO2

Fix the resonances at the abstract's rounded 9.24 and 4.70 eV. Weighted linear
least squares fits n^2-1 for the two oscillator strengths, using reported index
uncertainties propagated to n^2. This is OUR refit, not Uchida's original fitted
parameter table. Fits are intended for the original 400-1000 nm interval only.
"""
import json
import numpy as np

WAVELENGTH_UM = [.4047, .4358, .4678, .4800, .5086, .5461, .5893, .6328,
                 .6438, .690, .800, 1.00]
OBSERVED = {
    "o": [2.4315, 2.3834, 2.3478, 2.3366, 2.3150, 2.2931, 2.2738, 2.2597,
          2.2562, 2.2450, 2.226, 2.208],
    "e": [2.6157, 2.5583, 2.5164, 2.5036, 2.4779, 2.4520, 2.4295, 2.4119,
          2.4086, 2.3955, 2.373, 2.352],
}
UNCERTAINTY = {
    "o": [.0009, .0007, .0007, .0006, .0006, .0006, .0006, .0006,
          .0006, .0006, .0015, .0015],
    "e": [.0010, .0008, .0008, .0008, .0007, .0007, .0007, .0007,
          .0007, .0007, .0015, .0015],
}


def fit():
    w = np.asarray(WAVELENGTH_UM)
    resonances = [1.239841984/9.24, 1.239841984/4.70]
    X = np.column_stack([w*w/(w*w-a*a) for a in resonances])
    result = {"resonance_wavelengths_um": resonances, "branches": {}}
    for pol in "oe":
        n = np.asarray(OBSERVED[pol])
        weights = 1/(2*n*np.asarray(UNCERTAINTY[pol]))
        coefficients = np.linalg.lstsq(X*weights[:, None], (n*n-1)*weights, rcond=None)[0]
        residual = np.sqrt(1+X@coefficients)-n
        result["branches"][pol] = dict(strengths=coefficients.tolist(),
                                      max_abs_index_residual=float(np.max(np.abs(residual))),
                                      index_residuals=residual.tolist())
    return result


if __name__ == "__main__":
    print(json.dumps(fit(), indent=2))
