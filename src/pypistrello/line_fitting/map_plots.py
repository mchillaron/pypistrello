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
from scipy.interpolate import griddata

def make_a_map(
    table_results_fitting,
    cube_data,
    flux_col="line_flux_trapz",
    cmap="plasma",
    vmin=None,
    vmax=None,
    title="Line flux map"
):
    """
    Create and display a 2D map of line flux from a results table.

    Parameters
    ----------
    table_results_fitting : astropy.table.Table
        Table with columns 'x', 'y' (FITS coords, 1-based) and flux.
    nx, ny : int
        Spatial dimensions of the cube.
    flux_col : str
        Name of the flux column.
    cmap : str
        Matplotlib colormap.
    vmin, vmax : float, optional
        Color scale limits.
    title : str
        Plot title.

    Returns
    -------
    image : ndarray
        2D array with line flux values.
    """

    # Initialize image with NaNs
    nw, ny, nx = cube_data.shape
    image = np.full((ny, nx), np.nan)

    # Fill image
    for row in table_results_fitting:
        x = int(row["x"]) - 1   # FITS → numpy
        y = int(row["y"]) - 1
        image[y, x] = row[flux_col]

    # Plot
    plt.figure(figsize=(7, 6))
    im = plt.imshow(
        image,
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax
    )
    plt.colorbar(im, label=flux_col)
    plt.xlabel("X (pixel)")
    plt.ylabel("Y (pixel)")
    plt.title(title)
    plt.tight_layout()
    plt.show()

    return image


def build_2d_map(x, y, values, interpolate=True, method="nearest"):
    """
    Build a 2D map from x, y, values.

    Parameters
    ----------
    interpolate : bool
        If True, interpolate missing pixels
    method : str
        Interpolation method for griddata
    """

    x = np.asarray(x).astype(int)
    y = np.asarray(y).astype(int)

    nx = x.max() + 1
    ny = y.max() + 1

    zi = np.full((ny, nx), np.nan)

    if interpolate:
        xi = np.arange(nx)
        yi = np.arange(ny)
        xi_grid, yi_grid = np.meshgrid(xi, yi)

        zi = griddata(
            (x, y),
            values,
            (xi_grid, yi_grid),
            method=method
        )
    else:
        zi[y, x] = values

    return zi





def build_2d_map_always_interpolate(x, y, values):
    xi = np.unique(x)
    yi = np.unique(y)
    xi_grid, yi_grid = np.meshgrid(xi, yi)

    zi = griddata(
        (x, y),
        values,
        (xi_grid, yi_grid),
        method="nearest"
    )

    return xi_grid, yi_grid, zi

def save_contours(output_path, contours):
    np.savez(output_path, **contours)

def load_contours(contour_file):
    return dict(np.load(contour_file))
