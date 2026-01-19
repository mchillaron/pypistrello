#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Load an Astropy Table from FITS filepath."""

from astropy.table import Table

def load_fits_table(fits_path):
    """Load an Astropy Table from a FITS file.

    Parameters
    ----------
    fits_path : Path
        Path to the FITS file.

    Returns
    -------
    Table
        Astropy Table loaded from the FITS file.
    """

    table = Table.read(fits_path, format='fits')
    
    return table