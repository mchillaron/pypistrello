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

def run_powerbin(table_results_fitting, config_parameters):
    """This function uses the package Powerbin: 
    'PowerBin method by Cappellari (2025, MNRAS, 544, 1432)'  
    https://ui.adsabs.harvard.edu/abs/2025MNRAS.544.1432C
    """
    x = table_results_fitting["x"] - 1
    y = table_results_fitting["y"] - 1
    xy = np.column_stack([x, y])
    signal = table_results_fitting["line_flux_trapz"]
    noise = table_results_fitting["noise"]

    target_sn = config_parameters["target_sn"]
    additive = config_parameters.get("additive", True)

    if additive:
        # ADDITIVE CASE: (S/N)^2 is additive when noise is Poissonian.
        capacity_spec = (signal / noise)**2
    else:
        # NON-ADDITIVE CASE: Define a function for custom capacity logic.
        def capacity_spec(index):
            # Standard S/N for the bin
            sn = np.sum(signal[index]) / np.sqrt(np.sum(noise[index]**2))
            # Example of modelling correlated noise (commented out):
            # sn /= 1 + 1.07 * np.log10(len(index))
            return sn**2

    # Perform the binning. The target is target_sn**2 to match the capacity definition.
    pow = PowerBin(xy, capacity_spec, target_capacity=target_sn**2, verbose=1)

    # The binning is performed on (S/N)^2, but for plotting we use S/N.
    #pow.plot(capacity_scale='sqrt', ylabel='S/N')
    #plt.show(block=True)

    bin_num = pow.bin_num                           # its shape is (N,) where N is the number of spaxels in the original table, and each value is the bin ID assigned to that spaxel.
    bin_capacity = pow.bin_capacity                 # its shape is (M,) where M is the number of bins, and each value is the total capacity (S/N)^2 of that bin.
    print("There is a total number of bins of", np.max(bin_num))

    table_results_fitting["bin_id"] = bin_num
    table_results_fitting["bin_capacity"] = bin_capacity[bin_num]
    table_results_fitting["bin_snr"] = np.sqrt(bin_capacity[bin_num])
    print("Columns 'bin_id', 'bin_capacity' and 'bin_snr' have been added to the table with the results of the Voronoi binning")
    print("End of PowerBin Voronoi binning")

    return pow, table_results_fitting