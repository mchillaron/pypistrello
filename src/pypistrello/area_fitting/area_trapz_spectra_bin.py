#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#


from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

import numpy as np
import matplotlib.pyplot as plt

from ..file_loading.load_yaml_file import validate_region_config
from ..file_loading.save_table_fits import save_table_with_wcs_extension
from .window_regions import get_region_mask, apply_excluded_regions
from .continuum import fit_continuum
from .trapz import compute_line_flux
from .plotting import plot_trapz_spectrum
from .save_trapz_plots_pdf import save_trapz_plots_to_pdf
from .signal_to_noise import signal_to_noise


def area_trapz_spectra_bin(spectra, wavelength_range, config_parameters, 
                        analysis_table, redshift, line_restframe, debug_level):
    
    # Preparing parameters for the fitting process:
    line_obs = np.array(line_restframe) * (1+redshift)
    wavelength = wavelength_range
    zoom_limits = np.array(config_parameters["zoom_plot"])
    y_pad = config_parameters["y_padding"]
    
    poly_order_cont = config_parameters["poly_order_cont"]
    print(f"The continuum will be fitted with a polynomial of order {poly_order_cont}")

    validate_region_config(config_parameters)
    
    results_area = []
    total_spectra = spectra.shape[1]
    print(f"INFO: Area calculation will start for a total number of {total_spectra} spectra")
    
    # for loop reading spectrum by spectrum
    # dimensions of spectra are (n_lambda, n_bins)
    for bin in tqdm(range(spectra.shape[1]), desc="Integrating spectra (bins)", unit="bin"):

        print("Processing bin", bin+1, "out of", total_spectra)
        x_bin = analysis_table["x"][bin]
        y_bin = analysis_table["y"][bin]
        print("The coordinates of the bin are (x, y) =", (x_bin, y_bin))

        spectrum = spectra[:, bin]
        flux = spectrum          # extract one spectrum
        print(len(spectrum))  # debería dar 7341

        if np.all(flux == 0):    # Skip empty spectra if needed
            print("INFO: Skipping empty spectrum in bin", {bin+1})
            continue

        # Create a mask for the continuum:
        cont_mask, left_cont, right_cont = get_region_mask(
                wavelength,
                center=line_obs,
                window=config_parameters["window_continuum"],
                region=np.array(config_parameters["reg_continuum"]))
        
        # Prepare the regions to be excluded from cont_mask
        excluded_regions = []
        excl = config_parameters["reg_excluded"]
        if excl is not None:
            excluded_regions.extend(excl)
        
        # Line fit region must also be excluded from continuum
        fit_region = config_parameters["reg_fitting"]
        if fit_region is not None:
            excluded_regions.append(fit_region)

        if excluded_regions is not None:
            cont_mask = apply_excluded_regions(cont_mask, wavelength, excluded_regions)

        # Fitting the continuum
        cont_fit_func, coeffs, lambda_cont, flux_cont = fit_continuum(
            wavelength, flux, cont_mask,
            config_parameters["poly_order_cont"])

        # Line mask
        line_mask, line_left, line_right = get_region_mask(
            wavelength, center=line_obs,
            window=config_parameters["window_fitting"],
            region=config_parameters["reg_fitting"])

        area_trapz, lambda_line, flux_line_without_cont = compute_line_flux(
            wavelength, flux, line_mask, cont_fit_func,)

        # SNR (signal-to-noise ratio)
        line_snr, noise = signal_to_noise(wavelength, flux,
                                cont_mask, cont_fit_func,
                                line_mask, area_trapz)


        results_area.append((x_bin, y_bin, area_trapz, coeffs, noise, line_snr))

    results_area = np.array(results_area, dtype=object)

    n = len(analysis_table)
    analysis_table["bin_area_trapz"] = np.zeros(n)
    analysis_table["bin_cont_noise"] = np.zeros(n)
    analysis_table["bin_snr_trapz"] = np.zeros(n)
    analysis_table["bin_cont_coeffs"] = np.zeros((n, poly_order_cont + 1))

    # Add results to table making sure to match the order of the bins:
    for i, bin in enumerate(analysis_table["bin_id"]):
        mask = (
            (results_area[:, 0] == analysis_table["x"][i]) &
            (results_area[:, 1] == analysis_table["y"][i])
        )
        if np.sum(mask) == 1:
            area_trapz = results_area[mask, 2][0]
            coeffs = results_area[mask, 3][0] # This is an array of polynomial coefficients, not a numeric column
            noise = results_area[mask, 4][0]
            line_snr = results_area[mask, 5][0]

            analysis_table["bin_area_trapz"][i] = area_trapz
            analysis_table["bin_cont_noise"][i] = noise
            analysis_table["bin_snr_trapz"][i] = line_snr
            analysis_table["bin_cont_coeffs"][i] = coeffs
            
            
        else:
            print(f"WARNING: No unique match found for bin {bin} in results_area. Skipping this bin.")

    print("INFO: Area calculation completed and added to the table.")
    print(analysis_table[:5])
    return analysis_table
