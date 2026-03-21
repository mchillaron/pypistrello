#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np
from .fit_continuum_model import fit_continuum
from .line_models import gaussian_lmfit
from .line_models import fit_model_lmfit

def resolve_guess(value, default):
    """
    Resolve initial guess value.

    Parameters
    ----------
    value : user-provided value (can be None or 'auto')
    default : automatic value computed from data

    Returns
    -------
    float
    """
    if value is None or value == "auto":
        return default
    return value

def get_model_and_initial_params(config, x, y):
    # read the model to fit from the YAML
    model_name = config.get("model", "gaussian").lower()

    if model_name == "gaussian":
        model_func = gaussian_lmfit

        guesses = config.get("initial_guesses", {})
        p0 = {
            "amp": resolve_guess(guesses.get("amp"), np.nanmax(y)),
            "center": resolve_guess(guesses.get("center"), x[np.nanargmax(y)]),
            "sigma": resolve_guess(guesses.get("sigma"), 1.0),
        }

    else:
        raise ValueError(f"Model '{model_name}' not implemented")

    return model_func, p0

def fit_gaussian_spectrum_lmfit(wavelength, spectrum, config, debug=False):
    """
    Fit Gaussian to one spectrum using lmfit, with proper continuum handling.
    """

    # Fit continuum
    continuum_model = fit_continuum(wavelength, spectrum, config)

    if continuum_model is None:
        print("INFO: Continuum fit failed")
        return None

    continuum = continuum_model(wavelength)
    spectrum_sub = spectrum - continuum

    # Select fitting region
    
    lmin, lmax = config["reg_fitting"]
    mask = (wavelength >= lmin) & (wavelength <= lmax)

    x = wavelength[mask]
    y = spectrum_sub[mask]

    if len(x) < 5 or np.all(~np.isfinite(y)):
        raise ValueError("Not enough valid data points for fitting the line in 'reg_fitting'")

    model_func, p0 = get_model_and_initial_params(config, x, y)
    result = fit_model_lmfit(x, y, model_func, p0)

    if result is None:
        return None

    # Extract results
    amp = result.params["amp"].value
    mu = result.params["center"].value
    sigma = result.params["sigma"].value

    area = amp * sigma * np.sqrt(2 * np.pi)
    chi2 = result.chisqr
    residuals = result.residual

    # Debug plot
    if debug:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(wavelength, spectrum, label="Original")
        plt.plot(wavelength, continuum, label="Continuum")
        plt.plot(x, y, label="Line (cont sub)")
        plt.plot(x, result.best_fit, label="Gaussian fit")
        plt.legend()
        plt.title("DEBUG FIT")
        plt.show()

    return {
        "amp": amp,
        "mu": mu,
        "sigma": sigma,
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
    print(f"INFO: Starting Gaussian fitting for {n_spec} spectra")

    amp_arr = np.full(n_spec, np.nan)
    mu_arr = np.full(n_spec, np.nan)
    sigma_arr = np.full(n_spec, np.nan)
    area_arr = np.full(n_spec, np.nan)
    chi2_arr = np.full(n_spec, np.nan)

    for i in range(n_spec):

        if i % 50 == 0:
            print(f"INFO: Fitting spectrum {i}/{n_spec}")

        spec = spectra[:, i]

        result = fit_gaussian_spectrum_lmfit(
            wavelength,
            spec,
            config,
            debug=False
        )

        if result is None:
            continue

        amp_arr[i] = result["amp"]
        mu_arr[i] = result["mu"]
        sigma_arr[i] = result["sigma"]
        area_arr[i] = result["area"]
        chi2_arr[i] = result["chi2"]

    # Save results
    analysis_table["amp_gauss"] = amp_arr
    analysis_table["mu_gauss"] = mu_arr
    analysis_table["sigma_gauss"] = sigma_arr
    analysis_table["area_gauss"] = area_arr
    analysis_table["chi2_gauss"] = chi2_arr

    print("INFO: Gaussian fitting completed")

    return analysis_table