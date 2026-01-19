#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Load wavelength range from a CSV file as ndarray."""

import numpy as np

def load_wavelength_range(wavelength_path):
    """Load wavelength range from a CSV file to an NumPy array.

    Parameters
    ----------
    wavelength_path : Path
        Path to the CSV file containing the wavelength range.

    Returns
    -------
    ndarray
        Numpy array of wavelengths loaded from the CSV file.
    """

    wavelength_range = np.loadtxt(wavelength_path, delimiter=',')
    return wavelength_range