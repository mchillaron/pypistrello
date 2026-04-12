#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def compute_sim_snr(table, area_column, simulation_results): #real_table
        
    real_flux = table[area_column] # (L, 1) is 1D
    area_col_index = table.colnames.index(area_column)  #find the column index of the area_column in the table

    print("Computing SNR for column '", area_column, "' using the simulated measurements from the column index", area_col_index)
    flux_sim = simulation_results[:, :, area_col_index]      # (N_sim, L, 1) is 2D array with the line fluxes from all simulations
                                                # row: simulation, column: spectrum value
    noise = np.std(flux_sim, axis=0)            # (L, 1) is 1D array with the noise level at each line position for every separate spectrum, computed as the std of the fluxes from all simulations
    
    snr = real_flux / noise

    return snr