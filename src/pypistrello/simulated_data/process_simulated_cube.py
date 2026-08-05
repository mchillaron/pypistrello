#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from ..file_loading.load_fits_cube import read_fits_cube
from ..analysis_spectral_lines import main_trapz_fitting
from ..voronoi_binning.extract_spectra_from_table import extract_spectra_from_table
from ..voronoi_binning.build_voronoi_table import build_voronoi_table
from ..voronoi_binning.run_powerbin import run_powerbin
from ..voronoi_binning.sum_spectra_voronoi import sum_spectra_voronoi
from ..voronoi_binning.propagate_bin_to_spaxel import propagate_bin_to_spaxel_table
from ..analysis_tools.measure_spectra_properties import measure_spectra_properties
from ..area_fitting.area_trapz_spectra_bin import area_trapz_spectra_bin

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
        trapz_table,
        pow,
        pow_valid_mask, 
        ):

    print(f"INFO: Reading cube {cube_path.name} with data extension '{data_extension}'")
    primary_header, data_header, cube_data, wcs_info = read_fits_cube(cube_path, data_extension)

    if not real_cube_measured:

        print("INFO: Calculating trapezoidal areas for all spectra in this simulated cube")

        results_sim_table  = main_trapz_fitting(
            None,  # no output dir needed
            cube_data,
            wcs_info,
            wavelength_range,
            config_parameters,
            None,  # no table path needed
            redshift,
            line_restframe
        )

    else:
        
        print("INFO: Loading trapezoidal measurements")
        if trapz_table is None:
            raise RuntimeError("Trapzoidal measurements were not provided.")

        table_spaxels = trapz_table.copy()

        # check if voronoi has been carried out in the original cube
        if pow is not None:
            print("INFO: Applying the original Voronoi binning to simulated cube")

            n_spaxels = len(table_spaxels)

            bin_id = np.full(n_spaxels, -1, dtype=int)
            bin_center_x = np.full(n_spaxels, np.nan)
            bin_center_y = np.full(n_spaxels, np.nan)
            bin_capacity = np.full(n_spaxels, np.nan)
            bin_snr = np.full(n_spaxels, np.nan)

            bin_id[pow_valid_mask] = pow.bin_num                                    # its shape is (N,) where N is the number of spaxels in the original table, and each value is the bin ID assigned to that spaxel.
            bin_center_x[pow_valid_mask] = pow.xybin[pow.bin_num][:,0]
            bin_center_y[pow_valid_mask] = pow.xybin[pow.bin_num][:,1]
            bin_capacity[pow_valid_mask] = pow.bin_capacity[pow.bin_num]            # its shape is (M,) where M is the number of bins, and each value is the total capacity (S/N)^2 of that bin.
            bin_snr[pow_valid_mask] = np.sqrt(pow.bin_capacity[pow.bin_num])

            table_spaxels["bin_id"] = bin_id
            table_spaxels["bin_center_x"] = bin_center_x
            table_spaxels["bin_center_y"] = bin_center_y
            table_spaxels["bin_capacity"] = bin_capacity
            table_spaxels["bin_snr"] = bin_snr
            
            #pow_sim, table_spaxels = run_powerbin(table_spaxels, config_parameters, snr_table=snr_table)
            cube_binned_sim, bin_map_sim, cube_voronoi_sim = sum_spectra_voronoi(
                cube_data, table_spaxels, output_dir_path=None
            )
            analysis_table_sim = build_voronoi_table(table_spaxels, pow) # using the same pow from the original cube
            spectra_sim = cube_binned_sim   # (n_lambda, n_bins)
            analysis_table_sim = area_trapz_spectra_bin(spectra_sim, wavelength_range, config_parameters, 
                                                analysis_table_sim, redshift, line_restframe)   
            
        else:
            print("INFO: No Voronoi binning to simulated cube")
            spectra_sim = extract_spectra_from_table(cube_data, table_spaxels)
            analysis_table_sim = table_spaxels
            #pow_sim = None

        analysis_table_sim = measure_spectra_properties(spectra_sim, wavelength_range, config_parameters,
                                                        analysis_table_sim, redshift, line_restframe, real_cube_measured=real_cube_measured)
        if pow is not None:
            columns_to_copy = ["n_pix", "bin_center_x", "bin_center_y", 
                               "bin_area_trapz","bin_cont_noise","bin_snr_trapz","bin_cont_coeffs",  #"bin_snr_simulated",
                               "velocity", "offsets",
                               "amp_gauss", "mu_gauss", "sigma_gauss", "fwhm", "cont_gauss", "area_gauss", "chi2_gauss",
                               "amp_ha", "amp_nii6548", "amp_nii6583", "area_ha", "area_nii6548", "area_nii6583", "area_total,"
                               "amp1", "amp2", "mu1", "mu2", "area1", "area2", "area_total"]
            
            results_sim_table = propagate_bin_to_spaxel_table(table_spaxels, analysis_table_sim, columns_to_copy)
        else: # case pow_sim=None
            results_sim_table = analysis_table_sim
        
        
    print(results_sim_table)

    # the result of this function is an Astropy Table so we convert it to array
    results_sim, col_names = table_to_array(results_sim_table)

    return results_sim, col_names