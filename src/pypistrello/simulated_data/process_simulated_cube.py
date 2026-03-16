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

import numpy as np

def table_to_array(table):
    columns = ["x", "y", "line_flux_trapz", "noise", "line_snr"]  
    arr = np.column_stack([table[col] for col in columns])
    return arr

def process_simulated_cube(
        cube_path,
        wavelength_range,
        data_extension,
        config_parameters,
        redshift,
        line_restframe):

    print(f"INFO: Reading cube {cube_path.name} with data extension '{data_extension}'")
    primary_header, data_header, cube_data, wcs_info = read_fits_cube(cube_path, data_extension)

    print("INFO: Calculating trapezoidal areas for all spectra in this simulated cube")
    results_areatrapz_sim_table = main_line_fitting(
        None,  # no output dir needed
        cube_data,
        wcs_info,
        wavelength_range,
        config_parameters,
        None,  # no table path needed
        redshift,
        line_restframe
    )

    print(results_areatrapz_sim_table)
    results_areatrapz_sim = table_to_array(results_areatrapz_sim_table)

    return results_areatrapz_sim