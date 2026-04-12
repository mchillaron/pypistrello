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
from .line_models import gaussian_lmfit, gaussian_area_fixed_lmfit
from .line_models import fit_model_lmfit

def estimate_sigma(x, y):
    """
    Estimate sigma from second moment of the line.
    """

    y = y - np.nanmin(y)

    if np.sum(y) <= 0:
        return 1.0

    mean = np.sum(x * y) / np.sum(y)
    var = np.sum(y * (x - mean)**2) / np.sum(y)

    return np.sqrt(var)

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

def get_model_and_initial_params(config, x, y, wavelength, index, analysis_table):
    # read the model to fit from the YAML
    model_name = config.get("model", "gaussian").lower()

    if model_name == "gaussian":
        model_func = gaussian_lmfit
        
        if index == 0: # Only prints for the first spectrum
            print("INFO: The chosen model to fit is SIMPLE GAUSSIAN") 
            print("INFO: Reading initial guesses from configuration YAML file")
        
        guesses = config.get("initial_guesses", {})

        # better estimation of line center
        lambda_obs = config["line_restframe"][0] * (1 + config["redshift"])
        if analysis_table is not None and "offsets" in analysis_table.colnames:
            offset_pix = analysis_table["offsets"][index]
            dlambda = np.mean(np.diff(wavelength))
            mu_auto = lambda_obs + offset_pix * dlambda
        else:
            mu_auto = x[np.nanargmax(y)]

        # better estimation of sigma
        sigma_auto = estimate_sigma(x, y)
        sigma0 = resolve_guess(guesses.get("sigma"), sigma_auto)

        p0 = {
            "amp": resolve_guess(guesses.get("amp"), np.nanmax(y)),
            #"center": resolve_guess(guesses.get("center"), x[np.nanargmax(y)]),
            "center": resolve_guess(guesses.get("center"), mu_auto),
            "sigma": resolve_guess(guesses.get("sigma"), sigma_auto),
            #"sigma": resolve_guess(guesses.get("sigma"), 1.0),
        }
        print("INFO: initial guess model created before fitting")

    elif model_name == "gaussian_area_fixed":
        model_func = gaussian_area_fixed_lmfit
        if index == 0: # Only prints for the first spectrum
            print("INFO: The chosen model to fit is SIMPLE GAUSSIAN WITH FIXED AREA")
            print("INFO: Reading initial guesses from configuration YAML file")
        
        guesses = config.get("initial_guesses", {})

        # better estimation of line center
        lambda_obs = config["line_restframe"][0] * (1 + config["redshift"])

        if analysis_table is not None and "offsets" in analysis_table.colnames:
            offset_pix = analysis_table["offsets"][index]
            dlambda = np.mean(np.diff(wavelength))
            mu_auto = lambda_obs + offset_pix * dlambda
        else:
            mu_auto = x[np.nanargmax(y)]

        # better estimation of sigma
        sigma_auto = estimate_sigma(x, y)

        # Here we used the previuosly measured area
        # If 'bin_area_trapz' is present take this column; if not, we take 'area_trapz'
        if "bin_area_trapz" in analysis_table.colnames:
            area_fixed = analysis_table["bin_area_trapz"][index]
        else:
            area_fixed = analysis_table["area_trapz"][index]

        p0 = {
            "center": resolve_guess(guesses.get("center"), mu_auto),
            "sigma": resolve_guess(guesses.get("sigma"), sigma_auto),
            "area": area_fixed,  # fijo
        }
    else:
        raise ValueError(f"Model '{model_name}' not implemented")

    return model_func, p0

def fit_gaussian_spectrum_lmfit(wavelength, spectrum, config, index=None, analysis_table=None, debug=False):
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

    # new protection against non-finite values in x and y
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if len(x) < 5 or np.all(~np.isfinite(y)):
        print("INFO: The lengths of arrays x and y are: ",x.shape, y.shape)
        print("Number of not finite points in y:", np.sum(~np.isfinite(y)))
        raise ValueError("Not enough valid data points for fitting the line in 'reg_fitting'")

    model_func, p0 = get_model_and_initial_params(config, x, y, wavelength, index, analysis_table)
    result = fit_model_lmfit(x, y, model_func, p0, config)
    if debug and result is not None:
        print(result.params["sigma"])

    if result is None:
        # fill the dictionary results with NaN values if the fit failed
        return {
            "amp": np.nan,
            "mu": np.nan,
            "sigma": np.nan,
            "area": np.nan,
            "chi2": np.nan,
            "residuals": np.full_like(y, np.nan)
        }
        #return None

    # Extract results
    sigma = result.params["sigma"].value
    mu = result.params["center"].value

    if "amp" in result.params:
        amp = result.params["amp"].value
        area = amp * sigma * np.sqrt(2 * np.pi)
    else:
        # obtain amp from area and sigma if using area-fixed model
        area = result.params["area"].value
        amp = area / (sigma * np.sqrt(2*np.pi))

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
    config,
    offsets,
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
        spec = spectra[:, i]

        result = fit_gaussian_spectrum_lmfit(
            wavelength,
            spec,
            config,
            index=i,
            analysis_table=analysis_table,
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
    analysis_table["fwhm"] = 2.355 * sigma_arr
    analysis_table["area_gauss"] = area_arr
    analysis_table["chi2_gauss"] = chi2_arr

    print("INFO: Gaussian fitting completed")

    return analysis_table