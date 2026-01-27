#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from numina.array.wavecalib.crosscorrelation import periodic_corr1d

def crosscorrelate_spectra(cube_data, wavelength, config_parameters):

    nw, ny, nx = cube_data.shape

    x_ref, y_ref = np.array(config_parameters["crosscorr"]["reference_spectrum"]) # FITS format
    x_ref -= 1 # change to python indices
    y_ref -= 1
    spectrum_ref = cube_data[:, y_ref, x_ref] # correct to make it python-index

    crosscorr_region = np.array(config_parameters["continuum"]["continuum_region"])
    print(crosscorr_region)

    lambda_min_crosscorr, lambda_max_crosscorr = crosscorr_region
    print(lambda_min_crosscorr) # wavelength left
    print(lambda_max_crosscorr) # wavelength right

    # need to change from wavelength to pixel to be able to extract from spectra
    iw_min = np.searchsorted(wavelength, lambda_min_crosscorr)
    iw_max = np.searchsorted(wavelength, lambda_max_crosscorr)
    iw_min = max(iw_min, 0)
    iw_max = min(iw_max, nw)

    spectrum_region_reference = spectrum_ref[iw_min:iw_max]

    offsets_pixel = []
    fpeaks_croscorr = []

    for y in tqdm(range(ny), desc="Crosscorrelating spectra (rows)", unit="row"):
        for x in range(nx):
            spectrum_to_compare = cube_data[:, y, x]   # extract one spectrum
            if np.all(spectrum_to_compare == 0):       # Skip empty spectra if needed
                continue

            offset_pixel, fpeak_croscorr = periodic_corr1d(
                sp_reference=spectrum_region_reference,
                sp_offset=spectrum_to_compare[iw_min:iw_max],
                remove_mean=True,
                frac_cosbell=0.10,
                zero_padding=50,
                debugplot=0   # 0: no plots; 12: one plot (crosscorrelation only); 22: more plots
            )

            offsets_pixel.append(offset_pixel)
            fpeaks_croscorr.append(fpeak_croscorr)


    offsets_pixel_array = np.array(offsets_pixel)
    fpeak_croscorr_array = np.array(fpeaks_croscorr)

    return offsets_pixel_array, fpeak_croscorr_array


def convert_offset_velocity(offsets_pixel_array, wavelength,
                            redshift, line_restframe):
    
    dlambda_dp = np.mean(np.diff(wavelength)) # Å/pixel o nm/pixel
    lambda_obs = np.array(line_restframe) * (1 + redshift)

    # Considering linear wavelength axis:
    delta_lambda = offsets_pixel_array * dlambda_dp
    c_kms = 299792.458
    velocity_array = c_kms * delta_lambda / lambda_obs
    print(velocity_array)

    return velocity_array



