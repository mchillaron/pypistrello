#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def build_map_from_bins(bin_map, table, value_column):
    """
    Build a 2D image from Voronoi bin results.

    Parameters
    ----------
    bin_map : ndarray (ny, nx)
        bin_map[y, x] = bin_id

    table : astropy.table.Table (N_bins,)
        Table with one row per bin

    value_column : str
        Column to map (e.g. 'velocity', 'flux')

    Returns
    -------
    image : ndarray (ny, nx)
    """

    ny, nx = bin_map.shape

    image = np.full((ny, nx), np.nan)

    values = table[value_column]   # (N_bins,)

    valid = bin_map >= 0

    image[valid] = values[bin_map[valid]]

    print(f"[DEBUG] Built Voronoi map with shape {image.shape}")

    return image