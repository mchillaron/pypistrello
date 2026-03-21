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

    # Apply exclusion regions
    excluded_regions = []
    if excl is not None:
        excluded_regions.extend(excl)

    if fit_region is not None:
        excluded_regions.append(fit_region)

    if excluded_regions is not None:
        for reg in excluded_regions:
            mask &= ~((wavelength >= reg[0]) & (wavelength <= reg[1]))

    x = wavelength[mask]
    y = spectrum[mask]

    if len(x) < 5 or np.all(~np.isfinite(y)):
        return None

    order = config.get("poly_order_cont", 1)

    try:
        coeffs = np.polyfit(x, y, order)
        continuum_model = np.poly1d(coeffs)
        return continuum_model
    except Exception as e:
        print(f"WARNING: Continuum fit failed: {e}")
        return None