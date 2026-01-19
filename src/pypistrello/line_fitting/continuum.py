#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

# continuum.py
import numpy as np

def fit_continuum(
    wavelength,
    flux,
    cont_mask,
    poly_order):
    
    lam_cont = wavelength[cont_mask]
    flux_cont = flux[cont_mask]

    if len(lam_cont) < poly_order + 1:
        raise ValueError("Not enough points to fit continuum.")

    coeffs = np.polyfit(lam_cont, flux_cont, poly_order)
    return np.poly1d(coeffs), lam_cont, flux_cont
