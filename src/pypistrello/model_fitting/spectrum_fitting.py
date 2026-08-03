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
from .line_models import gaussian_lmfit, gaussian_area_fixed_lmfit, double_gaussian_lmfit, triplet_HaNII_lmfit
from .line_models import fit_model_lmfit
from .array_dictionary import FIELDS, COLUMN_NAMES, EMPTY_RESULTS

def estimate_sigma(x, y):
    """
    Estimate sigma from second moment of the line.
    """

    y = y - np.nanmin(y)    # remove negative values

    if np.sum(y) <= 0:
        return 1.0

    mean = np.sum(x * y) / np.sum(y)                # flux mass centroid
    var = np.sum(y * (x - mean)**2) / np.sum(y)     # pondered variance

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
            "center": resolve_guess(guesses.get("center"), mu_auto),
            "sigma": resolve_guess(guesses.get("sigma"), sigma_auto),
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
            "area": area_fixed,  # fixed area from trapezoidal integration
        }

    elif model_name == "double_gaussian":

        model_func = double_gaussian_lmfit
        if index == 0:
            print("INFO: Double Gaussian model selected")

        guesses = config.get("initial_guesses", {})

        lambda1 = config["line_fitting_restframe"][0]
        lambda2 = config["line_fitting_restframe"][1]
        print("INFO: Rest-frame wavelengths for double Gaussian fitting:", lambda1, lambda2)

        delta_lambda = lambda2 - lambda1            #separation between the two lines in rest-frame

        lambda_obs = lambda1 * (1 + config["redshift"])
        #lambda_obs2 = lambda2 * (1 + config["redshift"])
        #midpoint = 0.5 * (lambda_obs + lambda_obs2)

        if analysis_table is not None and "offsets" in analysis_table.colnames:
            offset_pix = analysis_table["offsets"][index]
            dlambda = np.mean(np.diff(wavelength))
            mu_auto = lambda_obs + offset_pix*dlambda
        else:
            mu_auto = x[np.nanargmax(y)]

        #mask_sigma = np.abs(x-mu_auto) < 3

        mask_sigma = (
            (x > lambda_obs - 3.0) &
            (x < lambda_obs + 3.0)
        )

        if np.sum(mask_sigma) > 5:
            sigma_auto = estimate_sigma(
                x[mask_sigma],
                y[mask_sigma]
            )
        else:
            sigma_auto = estimate_sigma(x, y)

        sigma_auto = estimate_sigma(
            x[mask_sigma],
            y[mask_sigma]
        )

        #sigma_auto = estimate_sigma(x, y)

        if "bin_area_trapz" in analysis_table.colnames:
            area_fixed = analysis_table["bin_area_trapz"][index]
        else:
            area_fixed = analysis_table["area_trapz"][index]

        p0 = {
            "center": resolve_guess(guesses.get("center"), mu_auto),
            "sigma": resolve_guess(guesses.get("sigma"), sigma_auto),
            "area": area_fixed,
            "amp2": resolve_guess(guesses.get("amp2"),np.nanmax(y)),
            "delta_lambda": delta_lambda
        }

    elif model_name == "triplet_hanii":

        model_func = triplet_HaNII_lmfit

        if index == 0:
            print("INFO: The chosen model is Halpha+[NII] triplet")

        guesses = config.get("initial_guesses", {})

        # better estimation of line center
        lambda_obs = config["line_restframe"][0] * (1 + config["redshift"])

        if analysis_table is not None and "offsets" in analysis_table.colnames:
            offset_pix = analysis_table["offsets"][index]
            dlambda = np.mean(np.diff(wavelength))
            mu_auto = lambda_obs + offset_pix * dlambda
        else:
            mu_auto = x[np.nanargmax(y)]

        sigma_auto = estimate_sigma(x, y)

        # For Halpha, we consider the area fixed for the gaussian, 
        # as the value calculated with trapezoids
        if "bin_area_trapz" in analysis_table.colnames:
            area_fixed = analysis_table["bin_area_trapz"][index]
        else:
            area_fixed = analysis_table["area_trapz"][index]

        p0 = {
            "center": resolve_guess(guesses.get("center"),mu_auto),
            "sigma": resolve_guess(guesses.get("sigma"),sigma_auto),
            "area": area_fixed,
            "amp_nii6583": resolve_guess(guesses.get("amp"), np.nanmax(y)/40), #maximum divided by 40 to start fitting from a low amplitude for NII, helps convergence
        }

    else:
        raise ValueError(f"Model '{model_name}' not implemented")

    return model_func, p0

def extract_gaussian_result(result):

    sigma = result.params["sigma"].value
    mu = result.params["center"].value

    gauss_norm = sigma * np.sqrt(2*np.pi)

    if "amp" in result.params:
        amp = result.params["amp"].value
        area = amp * gauss_norm
    else:
        area = result.params["area"].value
        amp = area / gauss_norm

    return {
        "amp": amp,
        "mu": mu,
        "sigma": sigma,
        "area": area,
        "chi2": result.chisqr,
        "residuals": result.residual,
    }

def extract_triplet_result(result):

    sigma = result.params["sigma"].value
    mu = result.params["center"].value

    gauss_norm = sigma * np.sqrt(2*np.pi)

    area_ha = result.params["area"].value
    amp_ha = area_ha / gauss_norm

    amp_nii6583 = result.params["amp_nii6583"].value
    amp_nii6548 = amp_nii6583 / 3.0

    area_nii6548 = amp_nii6548 * gauss_norm
    area_nii6583 = amp_nii6583 * gauss_norm

    return {

        "amp_ha": amp_ha,
        "amp_nii6548": amp_nii6548,
        "amp_nii6583": amp_nii6583,

        "mu": mu,
        "sigma": sigma,

        "area_ha": area_ha,
        "area_nii6548": area_nii6548,
        "area_nii6583": area_nii6583,
        "area_total": area_ha + area_nii6548 + area_nii6583,

        "chi2": result.chisqr,
        "residuals": result.residual,
    }

def extract_double_gaussian_result(result):
    
    sigma = result.params["sigma"].value
    mu1 = result.params["center"].value
    delta = result.params["delta_lambda"].value
    mu2 = mu1 + delta

    area1 = result.params["area"].value
    amp1 = area1/(sigma*np.sqrt(2*np.pi))

    amp2 = result.params["amp2"].value
    area2 = amp2*sigma*np.sqrt(2*np.pi)

    return {
        "amp1": amp1,
        "amp2": amp2,
        "mu1": mu1,
        "mu2": mu2,
        "sigma": sigma,
        "area1": area1,
        "area2": area2,
        "area_total": area1 + area2,
        "chi2": result.chisqr,
        "residuals": result.residual,
    }

RESULT_EXTRACTORS = {

    "gaussian": extract_gaussian_result,

    "gaussian_area_fixed": extract_gaussian_result,

    "triplet_hanii": extract_triplet_result,

    "double_gaussian": extract_double_gaussian_result,

}

def fit_gaussian_spectrum_lmfit(wavelength, spectrum, config, index=None, analysis_table=None, debug=False):
    """
    Fit Gaussian to one spectrum using lmfit, with proper continuum handling.
    """
    model = config.get("model", "gaussian").lower()

    # Fit continuum
    continuum_model = fit_continuum(wavelength, spectrum, config)

    if continuum_model is None:
        print("INFO: Continuum fit failed")
        return None

    continuum = continuum_model(wavelength)
    spectrum_sub = spectrum - continuum

    # Select fitting region
    #lmin, lmax = config["reg_fitting"]
    model_name = config.get("model", "gaussian").lower()

    if model in ("triplet_hanii", "double_gaussian"):
        lmin, lmax = config.get("reg_fitting_gaussians",
                                config["reg_fitting"])
    else:
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
        output = {
            key: np.nan
            for key in EMPTY_RESULTS[model]
        }
        output["residuals"] = np.full_like(y, np.nan)
        return output


    # Debug plot
    if debug:
        import matplotlib.pyplot as plt

        plt.figure()

        plt.plot(wavelength, spectrum, label="Original")
        plt.plot(wavelength, continuum, label="Continuum")
        plt.plot(x, y, label="Continuum subtracted")
        plt.plot(x, result.best_fit, label="Best fit")

        plt.legend()
        plt.show()

    return RESULT_EXTRACTORS[model](result)



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
    model = config.get("model", "gaussian").lower()
    arrays = {
        name: np.full(n_spec, np.nan)
        for name in FIELDS[model]
    }
    print(arrays.keys())

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

        for key in arrays:
            arrays[key][i] = result[key]


    # Save results
    for key, column in COLUMN_NAMES[model].items():
        analysis_table[column] = arrays[key]

    analysis_table["fwhm"] = 2.355 * arrays["sigma"]

    print("INFO: Gaussian fitting completed")

    return analysis_table