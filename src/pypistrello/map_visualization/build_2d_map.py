#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from scipy.interpolate import griddata
import numpy as np

def build_2d_map(x, y, values, interpolate=True, method="nearest"):
    """
    Build a 2D map from scattered (x, y) coordinates and associated values.

    This function converts tabular pixel coordinates into a 2D image.
    Depending on the ``interpolate`` flag, the map is either:
    - interpolated onto a regular grid using ``scipy.interpolate.griddata``, or
    - filled directly at the given pixel locations.

    Parameters
    ----------
    x : array-like
        X pixel coordinates (FITS convention, 1-based indexing).
    y : array-like
        Y pixel coordinates (FITS convention, 1-based indexing).
    values : array-like
        Data values associated with each (x, y) coordinate.
    interpolate : bool, optional
        If True, interpolate the scattered data onto a regular grid.
        If False, assign values directly to their pixel positions.
        Default is True.
    method : str, optional
        Interpolation method passed to ``scipy.interpolate.griddata``.
        Common options are ``"nearest"``, ``"linear"``, and ``"cubic"``.
        Default is "nearest".

    Returns
    -------
    zi : numpy.ndarray
        2D array representing the reconstructed map. Pixels without
        assigned data are filled with NaN.
    """

    x = np.asarray(x).astype(int) - 1 # convert from FITS to Python indexes
    y = np.asarray(y).astype(int) - 1

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