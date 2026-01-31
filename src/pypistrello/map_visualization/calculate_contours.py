#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np
import os

from .build_2d_map import build_2d_map

def load_contours(contour_file):
    """
    Load contour data saved as a NumPy pickle file.

    Parameters
    ----------
    contour_file : str
        Path to the .npy file containing contour data.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'levels'   : 1D numpy array of contour levels
        - 'segments' : list of lists of (N, 2) arrays with contour segments
    """
    contours = np.load(contour_file, allow_pickle=True)

    if isinstance(contours, np.ndarray): # np.load returns a numpy.ndarray with dtype=object
        contours = contours.item()

    if not isinstance(contours, dict):
        raise ValueError("Invalid contour file format: expected a dict")

    if "levels" not in contours or "segments" not in contours:
        raise ValueError("Contour file must contain 'levels' and 'segments'")

    return contours

def save_contours(output_path, contours):
    """
    Save contour information to a NumPy binary file.

    The contour data typically includes contour levels and the
    corresponding line segments. Storing them allows contours
    to be reused or overplotted without recomputing them.

    Parameters
    ----------
    output_path : str
        Path to the output ``.npy`` file where the contour data
        will be saved.
    contours : dict
        Dictionary containing contour information. Expected keys are:
        - 'levels': array-like
            Contour levels.
        - 'segments': list
            List of contour line segments as returned by Matplotlib.
    """
    np.save(output_path, contours, allow_pickle=True)

def calculate_contours(params, table, working_dir, map_choice, ax):
    """
    Compute and plot contour levels from tabular data.

    This function interpolates scattered data onto a 2D grid,
    estimates noise properties in a robust way, computes contour
    levels, plots them on the provided Matplotlib axis, and
    saves the contour information to disk.

    Parameters
    ----------
    params : dict
        Dictionary of parameters controlling contour computation,
        read from a YAML configuration file. 
    table : astropy.table.Table
        Input table containing at least the columns 'x', 'y',
        and the data column specified in ``params["data_column"]``.
        Coordinates are assumed to be in FITS pixel format.
    working_dir : str
        Path to the working directory where the contour file
        will be saved.
    map_choice : str
        Identifier of the current map, used to name the contour file.
    ax : matplotlib.axes.Axes
        Matplotlib axis object where the contours will be drawn.
    """
    x = table["x"] # FITS format from the Table Attention!
    y = table["y"] # FITS format from the Table
    data = table[params["data_column"]]

    zi_contours = build_2d_map(
        x, y, data,
        interpolate=True,
        method="linear"
    )

    finite = zi_contours[np.isfinite(zi_contours)] # work only with finite values
    if finite.size < 10:
        print("WARNING: Not enough finite pixels for contours")
        return

    # estimation of the distribution of data and robust deviation
    p_low = params.get("contour_low_percentile", 20)
    p_high = params.get("contour_high_percentile", 99)
    low_val = np.percentile(finite, p_low)
    high_val = np.percentile(finite, p_high)

    sigma = np.std(finite[finite < low_val])

    if not np.isfinite(sigma) or sigma <= 0:
        print("WARNING: Invalid sigma, using percentile-based contours")
        sigma = None

    # mask the noise using the estimated sigma above
    if sigma is not None:
        min_level = params.get("contour_min_sigma", 3) * sigma
    else:
        min_level = low_val

    zi_masked = np.where(zi_contours > min_level, zi_contours, np.nan)

    n_levels = params.get("contour_levels", 7) # Define levels for the contours
    if sigma is not None and min_level > 0 and high_val > min_level:
        levels = np.logspace(
            np.log10(min_level),
            np.log10(high_val),
            n_levels
        )
    else:
        levels = np.linspace(
            min_level,
            high_val,
            n_levels
        )

    contour_set = ax.contour(zi_masked, levels=levels, 
                             colors=params.get("contours_color", "black"),
                             linewidths=params.get("contours_linewidth", 1.0))

    contour_file = os.path.join(working_dir, f"{map_choice}_contours.npy")

    save_contours(
        contour_file,
        {
            "levels": contour_set.levels,
            "segments": contour_set.allsegs
        }
    )


def add_contours_to_plot(params, contour_file_loaded, ax):
    """
    Load precomputed contours from disk and overplot them on a map.

    This function reads contour levels and line segments from a saved
    NumPy file and plots them on the provided Matplotlib axis. It also
    allows selecting a subset of contour levels to display.

    Parameters
    ----------
    params : dict
        Dictionary of plotting parameters, read from a YAML file.
    contour_file_loaded : str
        Path to the ``.npy`` file containing the saved contour data.
    ax : matplotlib.axes.Axes
        Matplotlib axis object where the contours will be drawn.
    """
    contour_color = params.get("contours_color", "black")
    contour_alpha = params.get("contours_alpha", 0.9)
    contour_lw = params.get("contours_linewidth", 1)

    contours = load_contours(contour_file_loaded)

    levels = contours.get("levels")
    allsegs = contours.get("segments")
    if levels is None or allsegs is None:
        raise ValueError("Invalid contour file format: expected 'levels' and 'segments'")

    # Here we can choose also how many levels will be drawn
    level_indices = params.get("flux_contour_levels", None) # Example for the YAML: flux_contour_levels: [0, 2, 4]
    if level_indices is None:
        level_indices = range(len(levels))

    for i in level_indices:
        if i >= len(allsegs):
            continue
        for seg in allsegs[i]:
            ax.plot(
                seg[:, 0],
                seg[:, 1],
                color=contour_color,
                alpha=contour_alpha,
                linewidth=contour_lw
            )

    