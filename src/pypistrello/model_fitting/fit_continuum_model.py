#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def fit_continuum(wavelength, spectrum, config):
    """
    Fit continuum using polynomial over selected regions.
    """

    lmin, lmax = config["reg_continuum"]
    mask = (wavelength >= lmin) & (wavelength <= lmax)

    excl = config["reg_excluded"]
    fit_region = config["reg_fitting"]
    fit_region_gaussians = config["reg_fitting_gaussians"]

    # Apply exclusion regions
    excluded_regions = []
    if excl is not None:
        excluded_regions.extend(excl)

    if fit_region is not None:
        excluded_regions.append(fit_region)

    if fit_region_gaussians is not None:
        excluded_regions.append(fit_region_gaussians)

    if excluded_regions is not None:
        for reg in excluded_regions:
            mask &= ~((wavelength >= reg[0]) & (wavelength <= reg[1]))

    x = wavelength[mask]
    y = spectrum[mask]

    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if len(x) < 5:
        return None

    order = config.get("poly_order_cont", 1)

    try:
        coeffs = np.polyfit(x, y, order)

        if not np.all(np.isfinite(coeffs)):
            raise ValueError("Invalid coefficients")

        continuum_model = np.poly1d(coeffs)

    except Exception:
        # fallback orden 0
        try:
            coeffs = np.polyfit(x, y, 0)
            continuum_model = np.poly1d(coeffs)
        except Exception:
            # fallback mediana
            median_value = np.nanmedian(y)
            continuum_model = lambda w: np.full_like(w, median_value)

    return continuum_model
