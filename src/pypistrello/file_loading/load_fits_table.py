#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Load an Astropy Table from FITS filepath."""

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

def load_fits_table(fits_path):
    """
    Load an Astropy Table and associated WCS (if present) from a FITS file.
    """

    wcs = None

    with fits.open(fits_path) as hdul:

        # --- Load table
        if "RESULTS" in hdul:
            table = Table(hdul["RESULTS"].data)
        else:
            table = Table(hdul[1].data)

        # --- Load WCS from extension
        if "WCS_MAP" in hdul:
            wcs_header = hdul["WCS_MAP"].header
            wcs_try = WCS(wcs_header)

            if wcs_try.has_celestial:
                wcs = wcs_try.celestial
                print("INFO: Celestial WCS loaded from WCS_MAP")
                print(wcs)
            else:
                print("WARNING: WCS_MAP present but has no celestial axes")

    return table, wcs
