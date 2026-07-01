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
    #lmin, lmax = config["reg_fitting"]
    model_name = config.get("model", "gaussian").lower()

    if model_name == "triplet_hanii":
        lmin, lmax = config.get("reg_fitting_gaussians", config["reg_fitting"])

    elif model_name=="double_gaussian":
        lmin,lmax=config.get("reg_fitting_gaussians", config["reg_fitting"])

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
        # fill the dictionary results with NaN values if the fit failed
        if model_name == "triplet_hanii":
            return {
                "amp_ha": np.nan,
                "amp_nii6548": np.nan,
                "amp_nii6583": np.nan,
                "mu": np.nan,
                "sigma": np.nan,
                "area_ha": np.nan,
                "area_nii6548": np.nan,
                "area_nii6583": np.nan,
                "area_total": np.nan,
                "chi2": np.nan,
                "residuals": np.full_like(y, np.nan)
            }
        
        elif model_name == "double_gaussian":
            return {
                "amp1": np.nan,
                "amp2": np.nan,
                "mu1": np.nan,
                "mu2": np.nan,
                "sigma": np.nan,
                "area1": np.nan,
                "area2": np.nan,
                "area_total": np.nan,
                "delta_lambda": np.nan,
                "chi2": np.nan,
                "residuals": np.full_like(y, np.nan)
            }
    
        else:
            return {
                "amp": np.nan,
                "mu": np.nan,
                "sigma": np.nan,
                "area": np.nan,
                "chi2": np.nan,
                "residuals": np.full_like(y, np.nan)
            }

    # Extract results
    if config["model"].lower() == "triplet_hanii":
        sigma = result.params["sigma"].value
        mu = result.params["center"].value

        area_ha = result.params["area"].value
        amp_ha = area_ha / (sigma * np.sqrt(2*np.pi))
        
        #amp_ha = result.params["amp_ha"].value
        amp_nii6583 = result.params["amp_nii6583"].value
        amp_nii6548 = amp_nii6583 / 3.0  # fixed ratio for [NII] lines

        #area_ha = amp_ha * sigma * np.sqrt(2*np.pi)
        area_nii6548 = amp_nii6548 * sigma * np.sqrt(2*np.pi)
        area_nii6583 = amp_nii6583 * sigma * np.sqrt(2*np.pi)
        area_total = area_ha + area_nii6548 + area_nii6583

    elif config["model"].lower() == "double_gaussian":
        sigma = result.params["sigma"].value
        mu1 = result.params["center"].value
        delta = result.params["delta_lambda"].value
        mu2 = mu1 + delta
        area1 = result.params["area"].value
        amp1 = area1/(sigma*np.sqrt(2*np.pi))
        amp2 = result.params["amp2"].value
        area2 = amp2*sigma*np.sqrt(2*np.pi)
        area_total = area1 + area2

    else:
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

    if config["model"].lower() == "triplet_hanii":
        return {
            "amp_ha": amp_ha,
            "amp_nii6548": amp_nii6548,
            "amp_nii6583": amp_nii6583,
            "mu": mu,
            "sigma": sigma,
            "area_ha": area_ha,
            "area_nii6548": area_nii6548,
            "area_nii6583": area_nii6583,
            "area_total": area_total,
            "chi2": chi2,
            "residuals": residuals
        } 
    elif config["model"].lower()=="double_gaussian":
        return {
            "amp1": amp1,
            "amp2": amp2,
            "mu1": mu1,
            "mu2": mu2,
            "sigma": sigma,
            "area1": area1,
            "area2": area2,
            "area_total": area_total,
            "chi2": chi2,
            "residuals": residuals
        }
    
    else:
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

    if config["model"].lower() == "triplet_hanii":
        amp_ha_arr = np.full(n_spec, np.nan)
        amp_nii6548_arr = np.full(n_spec, np.nan)
        amp_nii6583_arr = np.full(n_spec, np.nan)
        area_ha_arr = np.full(n_spec, np.nan)
        area_nii6548_arr = np.full(n_spec, np.nan)
        area_nii6583_arr = np.full(n_spec, np.nan)

    elif config["model"].lower() == "double_gaussian":
        amp1_arr = np.full(n_spec, np.nan)
        amp2_arr = np.full(n_spec, np.nan)
        mu2_arr = np.full(n_spec, np.nan)
        area1_arr = np.full(n_spec, np.nan)
        area2_arr = np.full(n_spec, np.nan)

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
        
        if config["model"].lower() == "triplet_hanii":
            amp_ha_arr[i] = result["amp_ha"]
            amp_nii6548_arr[i] = result["amp_nii6548"]
            amp_nii6583_arr[i] = result["amp_nii6583"]
            mu_arr[i] = result["mu"]
            sigma_arr[i] = result["sigma"]
            area_ha_arr[i] = result["area_ha"]
            area_nii6548_arr[i] = result["area_nii6548"]
            area_nii6583_arr[i] = result["area_nii6583"]
            area_arr[i] = result["area_total"]
            chi2_arr[i] = result["chi2"]
        
        elif config["model"].lower() == "double_gaussian":
            amp1_arr[i] = result["amp1"]
            amp2_arr[i] = result["amp2"]
            mu_arr[i] = result["mu1"]  # or mu2, they are related by delta_lambda
            mu2_arr[i] = result["mu2"]
            sigma_arr[i] = result["sigma"]
            area1_arr[i] = result["area1"]
            area2_arr[i] = result["area2"]
            area_arr[i] = result["area_total"]
            chi2_arr[i] = result["chi2"]

        else: # gaussian or gaussian_area_fixed
            amp_arr[i] = result["amp"]
            mu_arr[i] = result["mu"]
            sigma_arr[i] = result["sigma"]
            area_arr[i] = result["area"]
            chi2_arr[i] = result["chi2"]

    # Save results
    if config["model"].lower() == "triplet_hanii":
        analysis_table["amp_ha"] = amp_ha_arr
        analysis_table["amp_nii6548"] = amp_nii6548_arr
        analysis_table["amp_nii6583"] = amp_nii6583_arr
        analysis_table["mu_gauss"] = mu_arr
        analysis_table["sigma_gauss"] = sigma_arr
        analysis_table["fwhm"] = 2.355 * sigma_arr
        analysis_table["area_ha"] = area_ha_arr
        analysis_table["area_nii6548"] = area_nii6548_arr
        analysis_table["area_nii6583"] = area_nii6583_arr
        analysis_table["area_total"] = area_arr
        analysis_table["chi2_gauss"] = chi2_arr

    elif config["model"].lower() == "double_gaussian":
        analysis_table["amp1"] = amp1_arr
        analysis_table["amp2"] = amp2_arr
        analysis_table["mu1_gauss"] = mu_arr
        analysis_table["mu2_gauss"] = mu2_arr
        analysis_table["sigma_gauss"] = sigma_arr
        analysis_table["fwhm"] = 2.355 * sigma_arr
        analysis_table["area1"] = area1_arr
        analysis_table["area2"] = area2_arr
        analysis_table["area_total"] = area_arr
        analysis_table["chi2_gauss"] = chi2_arr

    else:
        analysis_table["amp_gauss"] = amp_arr
        analysis_table["mu_gauss"] = mu_arr
        analysis_table["sigma_gauss"] = sigma_arr
        analysis_table["fwhm"] = 2.355 * sigma_arr
        analysis_table["area_gauss"] = area_arr
        analysis_table["chi2_gauss"] = chi2_arr

    print("INFO: Gaussian fitting completed")

    return analysis_table