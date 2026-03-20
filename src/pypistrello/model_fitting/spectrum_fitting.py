#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np
from .line_models import gaussian_lmfit
from .line_models import fit_model_lmfit

def fit_gaussian_spectrum_lmfit(wavelength, spectrum, config):
    """
    Fit a linemodel to one spectrum using lmfit.

    Returns a dictionary with fit parameters and statistics.
    """
    # Read from YAML the type of fitting we will apply:
    model_to_fit = config.get("model_to_fit", "gaussian").lower()
    if model_to_fit != "gaussian":
        print(f"WARNING: model_to_fit '{model_to_fit}' not recognized, defaulting to 'gaussian'")
        model_to_fit = "gaussian"
    
    if model_to_fit == "gaussian":
        chosen_model_func = gaussian_lmfit

    # Fitting window
    lmin, lmax = config["reg_continuum"]
    mask = (wavelength >= lmin) & (wavelength <= lmax)
    x = wavelength[mask]
    y = spectrum[mask]

    if len(x) < 5 or np.all(~np.isfinite(y)):
        return None

    # Initial guesses
    amp0 = np.nanmax(y) - np.nanmedian(y)
    mu0 = x[np.nanargmax(y)]
    sigma0 = 1.0
    cont0 = np.nanmedian(y)
    p0 = {"amp": amp0, "center": mu0, "sigma": sigma0, "cont": cont0}

    result = fit_model_lmfit(x, y, chosen_model_func, p0)
    if result is None:
        return None

    # Extract parameters
    amp = result.params["amp"].value
    mu = result.params["center"].value
    sigma = result.params["sigma"].value
    cont = result.params["cont"].value

    # Area of Gaussian
    area = amp * sigma * np.sqrt(2 * np.pi)

    # Chi2 and residuals
    chi2 = result.chisqr
    residuals = result.residual  # array same length as x

    return {
        "amp": amp,
        "mu": mu,
        "sigma": sigma,
        "cont": cont,
        "area": area,
        "chi2": chi2,
        "residuals": residuals
    }



def fit_gaussians_to_all_spectra_lmfit(
    spectra,
    wavelength,
    analysis_table,
    config
):
    """
    Loop over all spectra and fit Gaussian using lmfit.
    """

    n_spec = spectra.shape[1]
    print(f"INFO: Fitting Gaussian models to {n_spec} spectra with lmfit")

    # Prepare arrays to store results
    amp_arr = np.full(n_spec, np.nan)
    mu_arr = np.full(n_spec, np.nan)
    sigma_arr = np.full(n_spec, np.nan)
    cont_arr = np.full(n_spec, np.nan)
    area_arr = np.full(n_spec, np.nan)
    chi2_arr = np.full(n_spec, np.nan)
    residuals_list = [None] * n_spec # optionally, residuals could be stored in a list

    for i in range(n_spec):
        spec = spectra[:, i]

        result = fit_gaussian_spectrum_lmfit(
            wavelength,
            spec,
            config
        )

        if result is None:
            continue

        amp_arr[i] = result["amp"]
        mu_arr[i] = result["mu"]
        sigma_arr[i] = result["sigma"]
        cont_arr[i] = result["cont"]
        area_arr[i] = result["area"]
        chi2_arr[i] = result["chi2"]
        residuals_list[i] = result["residuals"]

    # Save into table
    analysis_table["amp_gauss"] = amp_arr
    analysis_table["mu_gauss"] = mu_arr
    analysis_table["sigma_gauss"] = sigma_arr
    analysis_table["cont_gauss"] = cont_arr
    analysis_table["area_gauss"] = area_arr
    analysis_table["chi2_gauss"] = chi2_arr
    # optionally, save residuals as object column
    analysis_table["residuals_gauss"] = residuals_list

    return analysis_table