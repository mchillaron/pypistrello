#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from ..file_loading.load_fits_cube import read_fits_cube
from ..analysis_spectral_lines import main_line_fitting
from ..voronoi_binning.extract_spectra_from_table import extract_spectra_from_table
from ..voronoi_binning.build_voronoi_table import build_voronoi_table
from ..voronoi_binning.run_powerbin import run_powerbin
from ..voronoi_binning.sum_spectra_voronoi import sum_spectra_voronoi
from ..voronoi_binning.propagate_bin_to_spaxel import propagate_bin_to_spaxel_table
from ..analysis_tools.measure_spectra_properties import measure_spectra_properties

import numpy as np

def table_to_array(table):
    columns = table.colnames
    print(f"The columns in the result table are {columns}") 
    arr = np.column_stack([table[col] for col in columns])
    return arr, columns

def process_simulated_cube(
        cube_path,
        wavelength_range,
        data_extension,
        config_parameters,
        redshift,
        line_restframe,
        real_cube_measured,
        snr_table,
        pow):

    print(f"INFO: Reading cube {cube_path.name} with data extension '{data_extension}'")
    primary_header, data_header, cube_data, wcs_info = read_fits_cube(cube_path, data_extension)

    if real_cube_measured:
        print("INFO: Calculating trapezoidal areas for all spectra in this simulated cube")
        # This is an Astropy table with "x", "y", "line_flux_trapz", "noise", "line_snr"
        table_spaxels = main_line_fitting(
            None,  # no output dir needed
            cube_data,
            wcs_info,
            wavelength_range,
            config_parameters,
            None,  # no table path needed
            redshift,
            line_restframe
        )
        
        # check if voronoi has been carried out in the original cube
        if pow is not None:
            print("INFO: Applying Voronoi binning to simulated cube")
            pow_sim, table_spaxels = run_powerbin(table_spaxels, config_parameters, snr_table=snr_table)
            cube_binned_sim, bin_map_sim, cube_voronoi_sim = sum_spectra_voronoi(
                cube_data, table_spaxels, output_dir_path=None
            )
            analysis_table_sim = build_voronoi_table(table_spaxels, pow_sim)
            spectra_sim = cube_binned_sim   # (n_lambda, n_bins)

        else: # case pow=None
            print("INFO: No Voronoi binning to simulated cube")
            spectra_sim = extract_spectra_from_table(cube_data, table_spaxels)
            analysis_table_sim = table_spaxels
            pow_sim = None

        analysis_table_sim = measure_spectra_properties(spectra_sim, wavelength_range, config_parameters,
                                                    analysis_table_sim, redshift, line_restframe)
        if pow_sim is not None:
            columns_to_copy = ["velocity","offsets","amp_gauss", "mu_gauss", "sigma_gauss", "fwhm", "cont_gauss", "area_gauss", "chi2_gauss"]
            results_sim_table = propagate_bin_to_spaxel_table(table_spaxels, analysis_table_sim, columns_to_copy)
        else: # case pow_sim=None
            results_sim_table = analysis_table_sim
        

    else:
        print("INFO: Calculating trapezoidal areas for all spectra in this simulated cube")
        results_sim_table = main_line_fitting(
            None,  # no output dir needed
            cube_data,
            wcs_info,
            wavelength_range,
            config_parameters,
            None,  # no table path needed
            redshift,
            line_restframe
        )

    print(results_sim_table)
    # the result of this function is an Astropy Table so we convert it to array
    results_sim, col_names = table_to_array(results_sim_table)

    return results_sim, col_names