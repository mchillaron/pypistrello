#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def propagate_bin_to_spaxel_table(spaxel_table, bin_table, columns_to_copy):
    """
    Propagate bin-level measurements back to spaxel-level table.

    Parameters
    ----------
    spaxel_table : astropy.table.Table (N_spaxels)
        Original table with 'bin_id'

    bin_table : astropy.table.Table (N_bins)
        Table with one row per bin

    columns_to_copy : list of str
        Columns from bin_table to propagate

    Returns
    -------
    spaxel_table : updated table
    """

    bin_id_spaxel = spaxel_table["bin_id"]   # (N_spaxels,)
    bin_id_bin = bin_table["bin_id"]         # (N_bins,)

    # Build mapping: bin_id → index in bin_table
    bin_index_map = {bid: i for i, bid in enumerate(bin_id_bin)}

    print(f"INFO: Propagating columns: {columns_to_copy}")

    for col in columns_to_copy:
        print("Processing column:", col)
        if col not in bin_table.colnames:
            print(f"WARNING: Column '{col}' not found in bin_table, skipping")
            continue

        values_bin = bin_table[col]   # (N_bins,)

        # Allocate output column
        #values_spaxel = np.zeros(len(spaxel_table))
        if col == "bin_cont_coeffs":
            values_spaxel = np.vstack([
                values_bin[bin_index_map[b]] for b in bin_id_spaxel
            ])
        else:
            values_spaxel = np.array([
                values_bin[bin_index_map[b]] for b in bin_id_spaxel
            ])

        for i in range(len(spaxel_table)):
            b = bin_id_spaxel[i]
            idx = bin_index_map[b]
            values_spaxel[i] = values_bin[idx]

        # Add to table
        spaxel_table[col] = values_spaxel

        print(f"INFO: Column '{col}' added to spaxel table")

    return spaxel_table