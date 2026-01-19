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
from .window_regions import get_region_mask, apply_excluded_regions
from .continuum import fit_continuum
from .trapz import compute_line_flux
from .plotting import plot_trapz_spectrum
from .save_trapz_plots_pdf import save_trapz_plots_to_pdf


def save_table(table, filename):
    table.write(filename, overwrite=True)

def main_line_fitting(output_dir_path, spectra_table, wavelength, 
                      config_parameters, table_parameters_path, 
                      redshift, line_restframe):
    
    
    line_obs = np.array(line_restframe) * (1+redshift)
    zoom_limits = config_parameters["plotting"]["zoom_plot"]
    y_pad = config_parameters["plotting"]["y_padding"]
    poly_order_cont = config_parameters["continuum"]["poly_order_cont"]

    validate_region_config(config_parameters)
    validate_region_config(config_parameters)

    results = []
    plot_inputs = []
    total_spectra=len(spectra_table)
    print(total_spectra)

    for row in tqdm(spectra_table, total=total_spectra, 
                    desc="Integrating area of spectrum", unit="spectrum"):
        flux = row["spec"]
        coords = (row["x"], row["y"])

        # Create a mask for the continuum before fitting it
        cont_mask, cont_left, cont_right = get_region_mask(
            wavelength,
            center=line_obs,
            window=config_parameters["continuum"]["window_cont"],
            region=config_parameters["continuum"]["continuum_region"])

        excluded_region = config_parameters["continuum"].get("excluded_region")
        if excluded_region is not None:
            cont_mask = apply_excluded_regions(cont_mask, wavelength, excluded_region)

        # Fitting the continuum
        cont_fit_func, lambda_cont, flux_cont = fit_continuum(
            wavelength, flux, cont_mask,
            config_parameters["continuum"]["poly_order_cont"])

        # --- Line mask
        line_mask, line_left, line_right = get_region_mask(
            wavelength, center=line_obs,
            window=config_parameters["line"]["window_line"],
            region=config_parameters["line"]["fit_region"])

        line_flux, lambda_line_sel, flux_line_sel = compute_line_flux(
            wavelength, flux, line_mask, cont_fit_func,)

        results.append((row["x"], row["y"], line_flux))

        plot_inputs.append(
            dict(
                wavelength=wavelength,
                flux=row["spec"],
                lambda_line_sel=lambda_line_sel,
                flux_line_sel=flux_line_sel,
                lambda_cont=lambda_cont,
                flux_cont=flux_cont,
                cont_fit_func=cont_fit_func,
                line_left=line_left,
                line_right=line_right,
                line_obs=line_obs,
                excluded_region=config_parameters["continuum"]["excluded_region"],
                zoom_limits=config_parameters["plotting"]["zoom_plot"],
                poly_order_cont=config_parameters["continuum"]["poly_order_cont"],
                coords=(row["x"], row["y"]),
                y_pad=config_parameters["plotting"]["y_padding"],
            )
        )

    table_results = Table(rows=results,names=("x", "y", "line_flux"))
    save_table(table_results, table_parameters_path)

    if config_parameters["plotting"]["save_pdf"]:
        pdf_path = output_dir_path / config_parameters["plotting"]["pdf_name"]
        save_trapz_plots_to_pdf(
            total_spectra,
            plot_inputs,
            pdf_path,
            *config_parameters["plotting"]["plots_per_page"],
        )

    return table_results
        
        