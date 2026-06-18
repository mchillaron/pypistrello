#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from astropy.io import fits
from cmap import Colormap

import matplotlib.pyplot as plt
import numpy as np

def plot_bin_map(bin_map):
    plt.figure(figsize=(6, 5))
    plt.imshow(bin_map, origin="lower", cmap='terrain')
    plt.colorbar(label="Bin ID")
    plt.title("Voronoi Bin Map")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()


def sum_spectra_voronoi(cube_data, table_results_fitting, output_dir_path, debug_level=0):
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

    Returns
    -------
    cube_binned : ndarray
        Binned spectra with shape (n_lambda, n_bins).
    var_binned : ndarray or None
        Binned variances, if var_cube is provided.
    """
    
    n_lambda, ny, nx = cube_data.shape                  # cube_data shape: (n_lambda, ny, nx)
    print(f"INFO: Cube shape: n_lambda={n_lambda}, ny={ny}, nx={nx}")

    # bin_map will store, for each pixel (y,x), the bin ID
    # shape: (ny, nx)
    # default = -1 means "no bin assigned"
    bin_map = np.full((ny, nx), -1, dtype=int)

    # Table coordinates are FITS (1-based), convert to Python (0-based)
    x = table_results_fitting["x"] - 1      # shape: (N_spaxels,)
    y = table_results_fitting["y"] - 1      # shape: (N_spaxels,)

    # bin_id per spaxel
    bin_id = table_results_fitting["bin_id"] # shape: (N_spaxels,)
    n_bins = np.max(bin_id) + 1             # Total number of bins
    bin_map[y, x] = bin_id                  # the map contains the bin ID for each spaxel, 
                                            # and -1 for spaxels not in the table (if any)

    print(f"INFO: Number of bins: {n_bins}")
    print(f"INFO: Number of spaxels in table: {len(bin_id)}")

    # cube_binned will contain ONE spectrum per bin
    # shape: (n_lambda, n_bins)

    cube_binned = np.zeros((n_lambda, n_bins))

    # Sum spectra per bin 
    for i in range(n_bins):
        mask = (bin_map == i)           # 2D mask (i is the bin ID value, we save in the mask the coordinates of all the spectra belonging to the same bin)
        n_pix_bin = np.sum(mask)        # Count number of pixels in this bin
        if debug_level == 2:
            print(f"DEBUG: Bin {i}: {n_pix_bin} pixels")

        spectra = cube_data[:, mask]    # (n_lambda, N_pix_bin) Extract all the spectra at mask coordinates from the data_cube
        spectra_sum = np.sum(spectra, axis=1) 
        cube_binned[:, i] = spectra_sum / n_pix_bin  # add to cube_binned the total spectrum after sum (n_lambda,) divided by the number of pixels in the bin (n_pix_bin) to get the average spectrum per bin


    print("INFO: Reconstruction of a data cube with the same dimensions as the original with Voronoi bins spectra")
    # cube_voronoi: same shape as original cube
    # but each pixel contains its bin spectrum

    cube_voronoi = np.full((n_lambda, ny, nx), np.nan)
    for i in range(n_bins):
        mask = (bin_map == i)                               # selects all the pixels belonging to the same bin
        cube_voronoi[:, mask] = cube_binned[:, i][:, None]  # cube_binned[:, i] is a 1D spectrum, [:, None] converts it into (n_lambda, 1)
                                                            # NumPy repeats it automatically over all the pixels in the same bin.

    print("INFO: Voronoi binning spectra completed")

    if debug_level >=1:
        plot_bin_map(bin_map)

    print("INFO: Saving Voronoi bin map to FITS file")
    if output_dir_path is not None:
        fits.writeto(
            output_dir_path / "bin_map.fits",
            bin_map,
            overwrite=True
        )

    return cube_binned, bin_map, cube_voronoi 


        



