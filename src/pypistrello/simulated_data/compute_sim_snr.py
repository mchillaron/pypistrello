#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def compute_sim_snr(real_table, simulation_results):

    real_flux = real_table["line_flux_trapz"]   # (L, 1) is 1D
    flux_sim = simulation_results[:, :, 2]      # (N_sim, L, 1) is 2D array with the line fluxes from all simulations
                                                # row: simulation, column: spectrum value
    noise = np.std(flux_sim, axis=0)            # (L, 1) is 1D array with the noise level at each line position for every separate spectrum, computed as the std of the fluxes from all simulations

    snr = real_flux / noise

    return snr