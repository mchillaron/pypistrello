#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def sum_spectra_voronoi(cube_data, table_results_fitting, var_cube=None):
    """
    Combine spectra in a 3D FITS cube according to Voronoi / PowerBin bins.

    Parameters
    ----------
    cube_data : ndarray
        Data cube with shape (n_lambda, ny, nx).
    table_results_fitting : astropy.table.Table
        Table containing at least columns:
        - 'x', 'y' (FITS coordinates, 1-based)
        - 'bin_id' (Voronoi bin index per spaxel)
    var_cube : ndarray, optional
        Variance cube with same shape as cube_data.

    Returns
    -------
    cube_binned : ndarray
        Binned spectra with shape (n_lambda, n_bins).
    var_binned : ndarray or None
        Binned variances, if var_cube is provided.
    """
    n_lambda, ny, nx = cube_data.shape

    bin_map = np.full((ny, nx), -1, dtype=int)

    x = table_results_fitting["x"] - 1   # FITS → Python
    y = table_results_fitting["y"] - 1

    bin_id = table_results_fitting["bin_id"]
    n_bins = np.max(bin_id) + 1
    bin_map[y, x] = bin_id                  # the map contains the bin ID for each spaxel, 
                                            # and -1 for spaxels not in the table (if any)

    cube_binned = np.zeros((n_lambda, n_bins))
    var_binned = np.zeros((n_lambda, n_bins)) if var_cube is not None else None

    # --- Sum spectra per bin ---
    for i in range(n_bins):
        mask = (bin_map == i)           # 2D mask
        spectra = cube_data[:, mask]    # (n_lambda, N_pix_bin)
        cube_binned[:, i] = np.sum(spectra, axis=1)

        if var_cube is not None:
            var_binned[:, i] = np.sum(var_cube[:, mask], axis=1)

    return cube_binned #additionally, var_binned if needed


        



