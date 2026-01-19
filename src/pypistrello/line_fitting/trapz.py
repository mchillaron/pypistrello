#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def compute_line_flux(
    wavelength,
    flux,
    line_mask,
    cont_fit_func):
    
    lam_sel = wavelength[line_mask]
    flux_sel = flux[line_mask]

    cont_fit = cont_fit_func(lam_sel)
    flux_line = flux_sel - cont_fit

    flux_int = np.trapezoid(flux_line, lam_sel)

    return flux_int, lam_sel, flux_line
