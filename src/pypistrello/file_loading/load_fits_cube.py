#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import FITSFixedWarning

import warnings


GREEN   = "\033[92m"
RESET   = "\033[0m"

def read_fits_cube(fits_path, ext):

    if ext == 0:
        print(f"{GREEN}INFO:{RESET} Using extension 0 as data_header")
    else:
        print(f"{GREEN}INFO:{RESET} Using extension 0 as primary_header and extension {ext} as data_header")


    with fits.open(fits_path) as hdul:
        warnings.simplefilter("ignore", FITSFixedWarning)
        primary_header = hdul[0].header if len(hdul) > 0 else None

        try:
            data_header = hdul[ext].header
            data = hdul[ext].data
        except IndexError:
            raise IndexError(f"Extension {ext} does not exist. File contains {len(hdul)} HDUs.")
        
        wcs_info = {
            "wcs": None,
            #"header": None,
            "naxis": None,
            "ra_cent": None,
            "dec_cent": None
        }

        if primary_header is not None:
            try:
                w = WCS(primary_header)
                if w.has_celestial:
                    wcs_info["wcs"] = w.celestial
                    #wcs_info["header"] = primary_header
                    wcs_header = primary_header
            except Exception:
                pass

        if wcs_info["wcs"] is None:
            try:
                w = WCS(data_header)
                if w.has_celestial:
                    wcs_info["wcs"] = w.celestial
                    #wcs_info["header"] = data_header
                    wcs_header = data_header
            except Exception:
                pass

        if wcs_info["wcs"] is not None:
            #hdr = wcs_info["header"]
            hdr = wcs_header

            nx = hdr.get("NAXIS1")
            ny = hdr.get("NAXIS2")
            if nx is not None and ny is not None:
                wcs_info["naxis"] = (nx, ny)

            # Apuntado central si existe
            wcs_info["ra_cent"] = hdr.get("RA", hdr.get("CRVAL1"))
            wcs_info["dec_cent"] = hdr.get("DEC", hdr.get("CRVAL2"))

    print("The WCS information extracted from the FITS headers is:")
    print(wcs_info)
    
    return primary_header, data_header, data, wcs_info
