#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def get_keyword(key, *headers, default=None, required=False):
    """
    Search for a FITS keyword in multiple headers (in order).

    Parameters
    ----------
    key : str
        FITS keyword to search for.
    headers : astropy.io.fits.Header
        Headers to search, in priority order.
    default : optional
        Default value if key is not found.
    required : bool
        If True, raises KeyError if keyword is not found.

    Returns
    -------
    value
        Keyword value.
    """
    for hdr in headers:
        if hdr is not None and key in hdr:
            return hdr[key]

    if required:
        raise KeyError(f"Keyword '{key}' not found in any provided FITS header.")

    return default


def get_wavelength_axis(data_header, primary_header):
    """
    Generate the wavelength axis from FITS WCS information.
    """

    # Priority: data_header → primary_header
    headers = (data_header, primary_header)

    crval3 = get_keyword("CRVAL3", *headers, required=True)
    crpix3 = get_keyword("CRPIX3", *headers, required=True)
    naxis3 = get_keyword("NAXIS3", *headers, required=True)
    
    # CDELT can be in different forms
    cdelt3 = get_keyword("CD3_3", *headers)
    if cdelt3 is None:
        cdelt3 = get_keyword("CDELT3", *headers, required=True)

    # FITS convention: axis indices start at 1
    indices = np.arange(1, naxis3 + 1)
    wavelength = crval3 + (indices - crpix3) * cdelt3

    return wavelength