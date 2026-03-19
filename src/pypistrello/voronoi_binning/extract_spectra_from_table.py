#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#
import numpy as np


def extract_spectra_from_table(cube_data, table):
    """
    Extract spectra following the exact order of the table.

    Parameters
    ----------
    cube_data : ndarray (n_lambda, ny, nx)
    table : astropy.table.Table with columns 'x', 'y' (1-based)

    Returns
    -------
    spectra : ndarray (n_lambda, N)
        Spectra ordered exactly like the table rows.
    """

    # Convert FITS → Python indexing
    x = table["x"] - 1   # shape (N,)
    y = table["y"] - 1   # shape (N,)

    # Extract spectra in correct order
    spectra = cube_data[:, y, x]   # shape (n_lambda, N)

    print(f"INFO: Extracted {spectra.shape[1]} spectra from cube")

    return spectra

def check_alignment(cube_data, table, spectra, n_checks=10):
    x = table["x"] - 1
    y = table["y"] - 1

    indices = np.random.choice(len(x), n_checks, replace=False)

    for i in indices:
        spec_from_cube = cube_data[:, y[i], x[i]]
        spec_from_array = spectra[:, i]

        if not np.allclose(spec_from_cube, spec_from_array):
            print(f"[ERROR] Mismatch at index {i}")
            return False

    print("INFO: Spectra alignment check PASSED")
    return True