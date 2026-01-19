#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#


"""Calculate the integrated spectrum over a user-defined spatial region.

    This function sums all individual spectra contained within a rectangular region 
    of the spatial plane defined by the coordinates (x1, y1) and (x2, y2). 
    The integrated spectrum can then be plotted for diagnostic purposes before 
    performing detailed line analysis."""

import numpy as np

def calculate_integrated_spectrum(spectra_table, wavelength_range, diagnostic_spectra):
    """Calculate the integrated spectrum over a user-defined spatial region.

    This function sums all individual spectra contained within a rectangular region 
    of the spatial plane defined by the coordinates (x1, y1) and (x2, y2). 
    The integrated spectrum can then be plotted for diagnostic purposes before 
    performing detailed line analysis.

    Parameters
    ----------
    spectra_table : astropy.table.Table
        An Astropy Table containing the spectral data. Each row must have at least:
        - 'x': X-coordinate of the spatial pixel
        - 'y': Y-coordinate of the spatial pixel
        - 'spec': 1D NumPy array containing the spectrum at that pixel
    wavelength_range : numpy.ndarray
        1D array of wavelength values (same length as each spectrum in `spectra_table`).
        The output spectrum will have the same length as this array.
    diagnostic_spectra : tuple of float
        Four coordinates defining the spatial region to integrate:
        (x1, x2, y1, y2), where (x1, y1) is the lower-left corner and
        (x2, y2) is the upper-right corner of the rectangular selection.

    Returns
    -------
    integrated_spectrum : numpy.ndarray
        1D array containing the sum of all spectra within the selected region. 
        The shape matches `wavelength_range`. If no spectra are found within the 
        region, a RuntimeError is raised.
    """

    x1, x2, y1, y2 = diagnostic_spectra

    integrated_spectrum = np.zeros_like(wavelength_range, dtype=float)
    n_spectra = 0
    for row in spectra_table:
        x = row["x"]
        y = row["y"]

        if x1 <= x <= x2 and y1 <= y <= y2:
            spectrum = row["spec"]
            if spectrum.shape != wavelength_range.shape:
                raise ValueError("Spectrum length does not match wavelength range.")
            
            integrated_spectrum += spectrum
            n_spectra += 1

    if n_spectra == 0:
        raise RuntimeError("No spectra found in selected region.")

    print(f"INFO: Integrated {n_spectra} spectra.")
    return integrated_spectrum