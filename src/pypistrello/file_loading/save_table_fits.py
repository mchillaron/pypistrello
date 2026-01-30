#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from astropy.io import fits

def save_table(table, filename):
    """Save an Astropy Table to disk.
    Parameters
    ----------
    table : astropy.table.Table
        The table to be saved.
    filename : str or pathlib.Path
        Output file path. The format is inferred from the file extension
        (e.g. ``.fits``, ``.ecsv``).
    
    """
    table.write(filename, overwrite=True)


def save_table_with_wcs(table, output_path, wcs=None): 
    """ Saves an Astropy Table to FITS, optionally adding WCS to the table header. 

    Parameters
    ----------
    table : astropy.table.Table
        Table containing the results to be written.
    output_path : str or pathlib.Path
        Path to the output FITS file.
    wcs : astropy.wcs.WCS or None, optional
        A 2D celestial WCS (e.g. RA/Dec) to be stored in the Table header.
        If ``None``, no WCS is written.
    """ 
    
    table_hdu = fits.BinTableHDU(table)
    if wcs is not None:
        wcs_header = wcs.to_header(relax=True)
        table_hdu.header.update(wcs_header)

    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        table_hdu
    ])

    hdul.writeto(output_path, overwrite=True)