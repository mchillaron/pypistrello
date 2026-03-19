#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np
from astropy.table import Table

def build_voronoi_table(table_results_fitting, pow):
    """
    Collapse spaxel-based table into a bin-based table.

    Parameters
    ----------
    table_results_fitting : astropy.table.Table
        Table with one row per spaxel. Must contain:
        - 'x', 'y' (1-based coordinates)
        - 'bin_id'

    pow : PowerBin object
        Contains Voronoi binning results.

    Returns
    -------
    bin_table : astropy.table.Table
        Table with one row per bin.
    """

    bin_id = table_results_fitting["bin_id"]    # shape: (N_spaxels,)
    unique_bins = np.unique(bin_id)             # shape: (N_bins,)

    n_bins = len(unique_bins)

    print(f"INFO: Building bin table")
    print(f"INFO: Number of bins: {n_bins}")

    # Prepare lists to store results:
    bin_ids = []
    x_list = []
    y_list = []
    n_pix_list = []
    snr_list = []

    for b in unique_bins:

        mask = (bin_id == b)   # shape: (N_spaxels,) ; mask selects all spaxels belonging to bin b

        # number of pixels in this bin
        n_pix = np.sum(mask)

        # use PowerBin centroid
        # pow.xybin[b] → (x_center, y_center) in Python coords (0-based)
        x_center_pow = pow.xybin[b, 0] + 1   # convert back to FITS (1-based)
        y_center_pow = pow.xybin[b, 1] + 1

        # SNR of the bin
        # pow.bin_capacity = (S/N)^2
        snr_bin = np.sqrt(pow.bin_capacity[b])

        #print(f"[DEBUG] Bin {b}: n_pix={n_pix}, SNR={snr_bin:.2f}")

        bin_ids.append(b)
        x_list.append(x_center_pow)   
        y_list.append(y_center_pow)

        n_pix_list.append(n_pix)
        snr_list.append(snr_bin)

    # Build output table
    bin_table = Table()

    bin_table["bin_id"] = bin_ids             # shape: (N_bins,)
    bin_table["x"] = x_list                   # xy position of the generators
    bin_table["y"] = y_list
    bin_table["n_pix"] = n_pix_list           # number of spaxels per bin
    bin_table["bin_snr"] = snr_list           # SNR per bin

    print("INFO: Bin table created successfully")
    print(bin_table[:5])

    return bin_table