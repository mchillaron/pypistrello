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

def fit_continuum(wavelength, flux, cont_mask, poly_order):
    
    """
    Fit a polynomial continuum to a selected region of a spectrum.

    This function selects wavelength and flux values corresponding to
    continuum regions, fits a polynomial of a given order to those points,
    and returns the resulting continuum model as a callable function.

    Parameters
    ----------
    wavelength : array-like
        Array of wavelength values.
    flux : array-like
        Array of flux values corresponding to `wavelength`.
    cont_mask : array-like of bool
        Boolean mask selecting the wavelength regions used to fit the continuum.
    poly_order : int
        Order of the polynomial used to fit the continuum.

    Returns
    -------
    cont_fit_func : numpy.poly1d
        Polynomial function representing the fitted continuum.
    lam_cont : array-like
        Wavelength values used for the continuum fit.
    flux_cont : array-like
        Flux values used for the continuum fit.

    Raises
    ------
    ValueError
        If there are not enough data points to fit the polynomial of the
        requested order.
    """

    lam_cont = wavelength[cont_mask]
    flux_cont = flux[cont_mask]

    if len(lam_cont) < poly_order + 1:
        raise ValueError("Not enough points to fit continuum.")

    coeffs = np.polyfit(lam_cont, flux_cont, poly_order)
    #print(f"The coefficients are {coeffs}")

    cont_fit_func = np.poly1d(coeffs)
    #print(f"The function to fit the continuum is {cont_fit_func}")
    
    return cont_fit_func, coeffs, lam_cont, flux_cont
