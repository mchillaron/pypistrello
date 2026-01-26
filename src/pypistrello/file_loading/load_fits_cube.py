#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from astropy.io import fits

def read_fits_cube(fits_path, ext):
    """
    Reads a FITS cube and returns the Primary header and the data
    from the first extension.

    Parameters
    ----------
    fits_path : str
        Path to the FITS file.
    ext : int
        Extension number in the FITS cube where data is found.
    Returns
    -------
    header : astropy.io.fits.Header
        Header of the FITS file.
    data : numpy.ndarray
        Data array from the first extension.
    """

    with fits.open(fits_path) as hdul:
        if ext == 0:
            primary_header = None
            data_header = hdul[0].header
        else:
            primary_header = hdul[0].header
            try:
                data_header = hdul[ext].header
            except IndexError:
                raise IndexError(
                    f"Extension {ext} does not exist. "
                    f"File contains {len(hdul)} HDUs."
                )

        data = hdul[ext].data

    return primary_header, data_header, data
