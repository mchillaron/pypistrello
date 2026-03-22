#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#


import numpy as np
from lmfit import Model

def gaussian_lmfit(x, amp, center, sigma):
    """Single Gaussian (continuum already removed)."""
    return amp * np.exp(-0.5 * ((x - center)/sigma)**2)

def fit_model_lmfit(x, y, model_func, p0, config):
    """
    Generic lmfit fitting function.

    Parameters
    ----------
    x : array
        Wavelength array
    y : array
        Spectrum values
    model_func : callable
        Function to fit (e.g., gaussian_lmfit)
    p0 : dict
        Initial guesses for parameters, keys must match model_func argument names

    Returns
    -------
    result : lmfit.model.ModelResult
    """
    model = Model(model_func)
    params = model.make_params(**p0)

    bounds = config.get("parameter_bounds", {})
    vary_flags = config.get("parameter_vary", {})

    for name, par in params.items():
        if name in bounds:              # Set bounds
            pmin, pmax = bounds[name]
            if pmin is not None:
                par.min = pmin
            if pmax is not None:
                par.max = pmax

        if name in vary_flags:         # Set vary flag
            par.vary = vary_flags[name]
        else:
            par.vary = True  # default: vary

    try:
        result = model.fit(y, params, x=x)
        return result
    except Exception as e:
        print("Fit failed:", e)
        return None