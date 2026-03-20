#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def compute_line_flux(wavelength, flux, line_mask, cont_fit_func):
    """
    Compute the integrated flux of a spectral line after continuum subtraction.

    This function extracts the wavelength and flux values corresponding to
    a spectral line region, subtracts a fitted continuum, and integrates
    the continuum-subtracted flux over the line region using the trapezoidal rule.

    Parameters
    ----------
    wavelength : array-like
        Array of wavelength values.
    flux : array-like
        Array of flux values corresponding to `wavelength`.
    line_mask : array-like of bool
        Boolean mask selecting the wavelength range of the spectral line.
    cont_fit_func : callable
        Function that evaluates the fitted continuum at given wavelengths.
        It must accept an array of wavelengths and return an array of
        continuum flux values.

    Returns
    -------
    flux_int : float
        Integrated line flux after continuum subtraction.
    lam_line : array-like
        Wavelength values within the line region.
    flux_line_no_cont : array-like
        Continuum-subtracted flux values within the line region.
    """
    
    lam_line = wavelength[line_mask]
    flux_line = flux[line_mask]

    cont_fit = cont_fit_func(lam_line)
    flux_line_no_cont = flux_line - cont_fit

    # integration of flux-cont over the line-region wavelengths
    flux_int = np.trapezoid(flux_line_no_cont, lam_line) 

    return flux_int, lam_line, flux_line_no_cont
