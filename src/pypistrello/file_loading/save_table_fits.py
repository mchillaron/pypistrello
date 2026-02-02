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


from astropy.io import fits
import numpy as np

def save_table_with_wcs_extension(
    table,
    output_path,
    wcs_info=None,
    overwrite=True
):
    """
    Save an Astropy Table to FITS and optionally add a WCS extension
    describing the reconstructed 2D map.

    Parameters
    ----------
    table : astropy.table.Table
        Result table to be written.
    output_path : str or Path
        Output FITS file.
    wcs_info : dict or None
        Dictionary containing WCS and related metadata.
        Expected keys:
          - 'wcs'     : astropy.wcs.WCS
          - 'naxis'   : (nx, ny)
          - 'ra_cent' : float (deg)
          - 'dec_cent': float (deg)
    """

    # --- Table HDU
    table_hdu = fits.BinTableHDU(table, name="RESULTS")
    hdus = [fits.PrimaryHDU(), table_hdu]

    # --- Optional WCS extension
    if wcs_info is not None and wcs_info.get("wcs") is not None:

        wcs = wcs_info["wcs"]
        wcs_header = wcs.to_header(relax=True)

        # Add extra metadata
        if wcs_info.get("ra_cent") is not None:
            wcs_header["RA"] = wcs_info["ra_cent"]
        if wcs_info.get("dec_cent") is not None:
            wcs_header["DEC"] = wcs_info["dec_cent"]

        wcs_header["EXTNAME"] = "WCS_MAP"

        # Crear data con shape correcto
        if wcs_info.get("naxis") is not None:
            nx, ny = wcs_info["naxis"]
            data = np.zeros((ny, nx), dtype=np.float32)
        else:
            data = np.zeros((1, 1), dtype=np.float32)

        wcs_hdu = fits.ImageHDU(
            data=data,
            header=wcs_header,
            name="WCS_MAP"
        )

        hdus.append(wcs_hdu)

    fits.HDUList(hdus).writeto(output_path, overwrite=overwrite)
