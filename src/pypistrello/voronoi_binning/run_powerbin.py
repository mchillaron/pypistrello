#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from powerbin import PowerBin

import matplotlib.pyplot as plt
import numpy as np

def run_powerbin(table_results_fitting, config_parameters,
                 debug_level=0, snr_table=None):
    """
    Perform Voronoi binning using PowerBin.

    PowerBin method by Cappellari (2025, MNRAS, 544, 1432)
    https://ui.adsabs.harvard.edu/abs/2025MNRAS.544.1432C
    """

    x = table_results_fitting["x"] - 1
    y = table_results_fitting["y"] - 1
    xy = np.column_stack([x, y])

    signal = np.asarray(table_results_fitting["area_trapz"], dtype=float)

    # Compute the noise
    if snr_table is None:
        print("No S/N table provided.")
        print("Using area_trapz and continuum noise.")
        noise = np.asarray(table_results_fitting["cont_noise"], dtype=float)

    else:
        print("S/N table provided.")
        print("Using simulated S/N values.")

        signal_to_noise = np.asarray(snr_table, dtype=float)

        noise = np.full_like(signal, np.nan) # noise will be computed only for valid S/N values, that is, spaxels with finite and positive SNR

        valid_noise = (
            np.isfinite(signal_to_noise)
            & (signal_to_noise > 0)
        )

        noise[valid_noise] = (
            signal[valid_noise] /
            signal_to_noise[valid_noise]
        )

    # Select valid spaxels
    valid = (
        np.isfinite(signal)
        & np.isfinite(noise)
        & (signal > 0)
        & (noise > 0)
    )

    n_total = len(signal)
    n_valid = np.sum(valid)

    print(f"Valid spaxels : {n_valid}/{n_total}")
    print(f"Rejected      : {n_total-n_valid}")

    signal = signal[valid]
    noise = noise[valid]
    xy = xy[valid]

    if debug_level > 0:
        print("\nPowerBin input summary")
        print("----------------------")
        print(f"Total spaxels           : {len(signal)}")
        print(f"Non-finite signal       : {np.sum(~np.isfinite(signal))}")
        print(f"Signal <= 0             : {np.sum(signal <= 0)}")
        print(f"Non-finite noise        : {np.sum(~np.isfinite(noise))}")
        print(f"Noise <= 0              : {np.sum(noise <= 0)}")
        print(f"Valid spaxels           : {n_valid}")

    
    # Capacity definition
    target_sn = config_parameters["target_sn"]
    additive = config_parameters.get("additive", True)

    if additive:
        # ADDITIVE CASE: (S/N)^2 is additive when noise is Poissonian.
        print("Computing additive capacity.")

        signal_to_noise = signal / noise
        capacity_spec = signal_to_noise**2

    else:
        # NON-ADDITIVE CASE: (np.sum(signal) / np.sqrt(np.sum(noise**2)))^2 This is the standard formula for uncorrelated noise
        print("Computing non-additive capacity.")

        def capacity_spec(index):
            sn = (
                np.sum(signal[index])
                / np.sqrt(np.sum(noise[index]**2))
            )
            return sn**2

    # Run powerbin
    # Perform the binning. The target is target_sn**2 to match the capacity definition.
    print(f"Running PowerBin with target S/N = {target_sn}")
    pow = PowerBin(xy, capacity_spec, target_capacity=target_sn**2, verbose=1)
    
    if debug_level > 0:
        pow.plot(capacity_scale="sqrt", ylabel="S/N")
        plt.show(block=True)

    # Recover results with original table size
    n_spaxels = len(table_results_fitting)

    bin_id = np.full(n_spaxels, -1, dtype=int)
    bin_center_x = np.full(n_spaxels, np.nan)
    bin_center_y = np.full(n_spaxels, np.nan)
    bin_capacity = np.full(n_spaxels, np.nan)
    bin_snr = np.full(n_spaxels, np.nan)
    #bin_center = np.full((n_spaxels, 2), np.nan)

    bin_id[valid] = pow.bin_num
    bin_center_x[valid] = pow.xybin[pow.bin_num][:,0]
    bin_center_y[valid] = pow.xybin[pow.bin_num][:,1]
    bin_capacity[valid] = pow.bin_capacity[pow.bin_num]
    bin_snr[valid] = np.sqrt(pow.bin_capacity[pow.bin_num])
    #bin_center[valid] = pow.xybin[pow.bin_num]

    table_results_fitting["bin_id"] = bin_id
    table_results_fitting["bin_center_x"] = bin_center_x
    table_results_fitting["bin_center_y"] = bin_center_y
    table_results_fitting["bin_capacity"] = bin_capacity
    table_results_fitting["bin_snr"] = bin_snr
    #table_results_fitting["bin_center"] = bin_center

    print(f"Number of Voronoi bins: {pow.xybin.shape[0]}")

    print(
        "Columns 'bin_id', 'bin_center', "
        "'bin_capacity' and 'bin_snr' added."
    )

    print("End of PowerBin Voronoi binning")
    print("\nColumn lengths:")
    for name in table_results_fitting.colnames:
        print(name, len(table_results_fitting[name]))

    return pow, table_results_fitting, valid


def run_powerbin_noneg(table_results_fitting, config_parameters, debug_level=0, snr_table=None):
    """This function uses the package Powerbin: 
    'PowerBin method by Cappellari (2025, MNRAS, 544, 1432)'  
    https://ui.adsabs.harvard.edu/abs/2025MNRAS.544.1432C
    """
    x = table_results_fitting["x"] - 1
    y = table_results_fitting["y"] - 1
    xy = np.column_stack([x, y])

    if snr_table is None:
        print("No S/N table provided, using area_trapz and noise from the fitting results to compute the capacity.")
        signal = table_results_fitting["area_trapz"]
        noise = table_results_fitting["cont_noise"]

        signal_to_noise = signal / noise
        n_signal_neg = np.sum(signal_to_noise < 0)
        signal[signal < 0] = 0.0
        print(f"Spaxels with negative signal: {n_signal_neg}")
        
    else:
        print("S/N table provided, using area_trapz and S/N to compute the capacity.")
        signal = table_results_fitting["area_trapz"]
        signal_to_noise = snr_table # In this case we obtained the noise from simulations
        noise = signal / signal_to_noise

        n_signal_neg = np.sum(signal_to_noise < 0)
        signal[signal < 0] = 0.0
        print(f"Spaxels with negative signal: {n_signal_neg}")

    target_sn = config_parameters["target_sn"]
    additive = config_parameters.get("additive", True)

    if additive:
        print("Computing capacity as addittive")
        # ADDITIVE CASE: (S/N)^2 is additive when noise is Poissonian.
        signal_to_noise = signal / noise
        capacity_spec = signal_to_noise**2
    else:
        # NON-ADDITIVE CASE: (np.sum(signal) / np.sqrt(np.sum(noise**2)))^2 This is the standard formula for uncorrelated noise
        print("Computing capacity as non addittive")
        sn_values = []
        def capacity_spec(index):
            sn = np.sum(signal[index]) / np.sqrt(np.sum(noise[index]**2))
            sn_values.append(sn)
            return sn**2   

    # Perform the binning. The target is target_sn**2 to match the capacity definition.
    pow = PowerBin(xy, capacity_spec, target_capacity=target_sn**2, verbose=1)
    
    # The binning is performed on (S/N)^2, but for plotting we use S/N.
    if debug_level > 0:
        pow.plot(capacity_scale='sqrt', ylabel='S/N')
        plt.show(block=True)

    bin_num = pow.bin_num                           # its shape is (N,) where N is the number of spaxels in the original table, and each value is the bin ID assigned to that spaxel.
    bin_capacity = pow.bin_capacity                 # its shape is (M,) where M is the number of bins, and each value is the total capacity (S/N)^2 of that bin.
    bin_center = pow.xybin                          # shape is (M, 2) where M is the number of bins, and each row is the center coordinates of that bin.
    print("There is a total number of bins of", np.max(bin_num))

    table_results_fitting["bin_id"] = bin_num
    table_results_fitting["bin_center"] = bin_center[bin_num]
    table_results_fitting["bin_capacity"] = bin_capacity[bin_num]
    table_results_fitting["bin_snr"] = np.sqrt(bin_capacity[bin_num])
    
    print("Columns 'bin_id', 'bin_center', 'bin_capacity' and 'bin_snr' have been added to the table with the results of the Voronoi binning")
    print("End of PowerBin Voronoi binning")

    return pow, table_results_fitting