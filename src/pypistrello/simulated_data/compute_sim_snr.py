#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np
import matplotlib.pyplot as plt

def compute_sim_snr(table, area_column, simulation_results, config_parameters, debug_level=0): #real_table
        
    real_flux = table[area_column] # (L, 1) is 1D
    print("The real flux values for SNR are: ", real_flux)
    area_col_index = table.colnames.index(area_column)  # Find the column index of the area_column in the table
    print("The column {area_column} is at index {area_col_index} in the table".format(area_column=area_column, area_col_index=area_col_index))
    k_sigma = config_parameters.get("k_sigma", 1)       # Get k_sigma from config, default to 1 if not provided

    print("Computing SNR for column '", area_column, "' using the simulated measurements from the column index", area_col_index)
    flux_sim = simulation_results[:, :, area_col_index]     # (N_sim, L, 1) is 2D array with the line fluxes from all simulations
                                                            # row: simulation, column: spectrum value

    if k_sigma == 0:
        noise = np.std(flux_sim, axis=0)                    # (L, 1) is 1D array with the noise level at each line position for every separate spectrum, computed as the std of the fluxes from all simulations
    elif k_sigma == 1:
        p16 = np.percentile(flux_sim, 16, axis=0)
        p84 = np.percentile(flux_sim, 84, axis=0)
        noise = (p84 - p16) / 2     
    elif k_sigma ==2:
        p5 = np.percentile(flux_sim, 5, axis=0)
        p95 = np.percentile(flux_sim, 95, axis=0)
        noise = (p95 - p5) / 2
    elif k_sigma == 3:
        p1 = np.percentile(flux_sim, 1, axis=0)
        p99 = np.percentile(flux_sim, 99, axis=0)
        noise = (p99 - p1) / 2                                         

    snr = real_flux / noise

    if debug_level >0:
        plot_snr_debug(flux_sim, real_flux)

    return snr

def plot_snr_debug(flux_sim, real_flux):
    N_sim, L = flux_sim.shape
    n_plots=10
    random_indices = np.random.choice(L, size=min(n_plots, L), replace=False)
    for i, idx in enumerate(random_indices):
        plt.figure()
        
        values = flux_sim[:, idx]
        plt.hist(values, bins=30, color='C2', alpha=0.7, label='Simulated Fluxes')
        
        # Statistics
        std_val = np.std(values)
        mean_val = np.mean(values)
        real_val = real_flux[idx]
        p16 = np.percentile(values, 16)
        p84 = np.percentile(values, 84)
        p5 = np.percentile(values, 5)
        p95 = np.percentile(values, 95)
        
        plt.axvline(mean_val, linestyle='--', label=f"mean = {mean_val:.3e}")
        plt.axvline(mean_val + std_val, linestyle='-', color='mediumaquamarine', label=f"+1 std = {std_val:.3e}")
        plt.axvline(mean_val - std_val, linestyle='-', color='mediumaquamarine', label=f"-1 std")
        plt.axvspan(mean_val - std_val, mean_val + std_val, color='aquamarine',alpha=0.15, label="±1 std")

        plt.axvline(p16, linestyle='-', color='plum', label=f"p16 = {p16:.3e}")
        plt.axvline(p84, linestyle='-', color='plum', label=f"p84 = {p84:.3e}")

        plt.axvline(p5, linestyle='-', color='lightcoral', label=f"p5 = {p5:.3e}")
        plt.axvline(p95, linestyle='-', color='lightcoral', label=f"p95 = {p95:.3e}")

        plt.axvline(real_val, linestyle='-', linewidth=2.5, color='C1',label=f"real = {real_val:.3e}")
        
        plt.title(f"Simulated Flux distribution (index{idx})")
        plt.xlabel("Flux")
        plt.ylabel("Frequency")
        plt.legend()
        plt.show()