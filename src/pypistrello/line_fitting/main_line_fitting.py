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

def main_line_fitting(output_dir_path, cube_data, wcs_info,
                      wavelength, config_parameters, 
                      table_parameters_path, redshift, line_restframe):
    
    # Preparing parameters for the fitting process:
    line_obs = np.array(line_restframe) * (1+redshift)

    zoom_limits = np.array(config_parameters["zoom_plot"])
    y_pad = config_parameters["y_padding"]
    
    poly_order_cont = config_parameters["poly_order_cont"]
    print(f"The continuum will be fitted with a polynomial of order {poly_order_cont}")

    validate_region_config(config_parameters)
    
    results_area = []
    plot_inputs = []
    nw, ny, nx = cube_data.shape
    total_spectra = nx * ny
    print(f"INFO: Area calculation will start for a total number of {total_spectra} spectra")
    print("The steps are:")
    print("Extracting the continuum mask without excluding no regions")
    print(f"The exclusion regions will be substracted from continuum mask")
    print(f"Fitting the continuum using a polynomial of order {config_parameters["poly_order_cont"]}")
    print("Extracting the line mask from wavelength range")
    print("Integrating the area of the line with Trapezoids")
    print("Calculating SNR")

    for y in tqdm(range(ny), desc="Integrating spectra (rows)", unit="row"):
        for x in range(nx):

            flux = cube_data[:, y, x]   # extract one spectrum
            if np.all(flux == 0):       # Skip empty spectra if needed
                continue

            x_fits = x + 1              # Convert numpy indices → FITS coordinates (1-based)
            y_fits = y + 1
            coords = (x_fits, y_fits)

            # Create a mask for the continuum
            cont_mask, cont_left, cont_right = get_region_mask(
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
            cont_fit_func, lambda_cont, flux_cont = fit_continuum(
                wavelength, flux, cont_mask,
                config_parameters["poly_order_cont"])

            # Line mask
            line_mask, line_left, line_right = get_region_mask(
                wavelength, center=line_obs,
                window=config_parameters["window_fitting"],
                region=config_parameters["reg_fitting"])

            line_flux_trapz, lambda_line, flux_line_without_cont = compute_line_flux(
                wavelength, flux, line_mask, cont_fit_func,)

            # SNR (signal-to-noise ratio)
            line_snr, noise = signal_to_noise(wavelength, flux,
                                  cont_mask, cont_fit_func,
                                  line_mask, line_flux_trapz)


            results_area.append((x_fits, y_fits, line_flux_trapz, noise, line_snr))

            #the plotinput info was here

    # Convert list of tuples → columns
    results_area = np.array(results_area)

    x_fits = results_area[:, 0].astype(int)
    y_fits = results_area[:, 1].astype(int)
    line_flux_trapz = results_area[:, 2].astype(float)
    noise = results_area[:, 3].astype(float)
    line_snr = results_area[:, 4].astype(float)

    print(f"Writing the line AREA and SNR to FITS table")
    table_results_fitting = Table(
        [x_fits, y_fits, line_flux_trapz, noise, line_snr],
        names=("x", "y", "line_flux_trapz", "noise", "line_snr")
    )
    table_results_fitting.meta["LINE"] = config_parameters["line_name"]
    #table_results_fitting.meta["LINE_CEN"] = line_obs. #Attribute `LINE_CEN` of type <class 'numpy.ndarray'> cannot be added to FITS Header
    #table_results_fitting.meta["FLUX_UNIT"] = "erg/s/cm2/AA" #VerifyWarning: Keyword name 'FLUX_UNIT' is greater than 8 characters or contains characters not allowed by the FITS standard;
    table_results_fitting.meta["COMMENT"] = "Integrated line flux after continuum subtraction"
    table_results_fitting.meta["X_AXIS"] = "pixel axis 1"
    table_results_fitting.meta["Y_AXIS"] = "pixel axis 2"

    if table_parameters_path is not None:
        print(f"Saving FITS table to {table_parameters_path}")
        save_table_with_wcs_extension(
            table_results_fitting,
            table_parameters_path,
            wcs_info=wcs_info
        )

    if config_parameters["save_pdf"]:
        
        print("Creting a PDF file with analysis plots")
        pdf_path = output_dir_path / config_parameters["pdf_name"]

        print("Collapsing information from all spectra in this cube to prepare the plots")
        plot_inputs.append(
                dict(
                    wavelength=wavelength,
                    flux=flux,
                    lambda_line=lambda_line,
                    flux_line=flux[line_mask],   #the flux in the line region, including continuum level
                    lambda_cont=lambda_cont,
                    flux_cont=flux_cont,
                    cont_fit_func=cont_fit_func,
                    line_left=line_left,
                    line_right=line_right,
                    line_obs=line_obs,
                    excluded_region=config_parameters["reg_excluded"],
                    zoom_limits=config_parameters["zoom_plot"],
                    poly_order_cont=config_parameters["poly_order_cont"],
                    coords=coords,
                    y_pad=config_parameters["y_padding"],
                )
            )
        
        save_trapz_plots_to_pdf(
            total_spectra,
            plot_inputs,
            pdf_path,
            *config_parameters["plots_per_page"],
        )

    return table_results_fitting
        
        