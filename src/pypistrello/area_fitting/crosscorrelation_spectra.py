#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np
from tqdm import tqdm

from numina.array.wavecalib.crosscorrelation import periodic_corr1d

def plot_offsets(offsets):

    import matplotlib.pyplot as plt

    plt.figure()
    plt.hist(offsets[~np.isnan(offsets)], bins=100)
    plt.xlabel("Pixel offset")
    plt.ylabel("N")
    plt.title("Offset distribution")
    plt.show()

def crosscorrelate_spectra_unified(
        spectra,
        wavelength,
        config_parameters,
        analysis_table,
        debug_level):
    """
    Cross-correlate all spectra against a reference spectrum.

    Parameters
    ----------
    spectra : ndarray (n_lambda, N)
        Spectra (spaxels OR Voronoi bins)

    wavelength : ndarray (n_lambda,)
        Wavelength array

    analysis_table : astropy.table.Table (length N)
        Must contain 'x', 'y' (1-based coordinates)

    Returns
    -------
    offsets : ndarray (N,)
    fpeaks  : ndarray (N,)
    """

    import numpy as np
    from tqdm import tqdm

    n_lambda, N = spectra.shape

    print(f"INFO: Cross-correlation on {N} spectra")

    #  Get reference spectrum from coordinates
    x_ref, y_ref = np.array(config_parameters["ref_spec"])  # FITS coords

    # Find closest object in the table
    dx = analysis_table["x"] - x_ref
    dy = analysis_table["y"] - y_ref

    dist = np.sqrt(dx**2 + dy**2)
    ref_index = np.argmin(dist)

    print(f"INFO: Reference spectrum index: {ref_index}")
    print(f"INFO: Closest position: x={analysis_table['x'][ref_index]}, y={analysis_table['y'][ref_index]}")

    spectrum_ref = spectra[:, ref_index]

    # Define wavelength region
    
    lambda_min, lambda_max = config_parameters["reg_continuum"]

    iw_min = np.searchsorted(wavelength, lambda_min)
    iw_max = np.searchsorted(wavelength, lambda_max)

    spectrum_region_ref = spectrum_ref[iw_min:iw_max]

    print(f"[DEBUG] Wavelength window: {lambda_min} - {lambda_max}")
    print(f"[DEBUG] Pixel range: {iw_min}:{iw_max}")

    # Loop over spectra
    offsets = np.zeros(N)
    fpeaks = np.zeros(N)

    for i in tqdm(range(N), desc="Cross-correlating"):

        spec = spectra[:, i]

        # skip bad spectra
        if np.all(~np.isfinite(spec)):
            offsets[i] = np.nan
            fpeaks[i] = np.nan
            continue

        off, fp = periodic_corr1d(
            sp_reference=spectrum_region_ref,
            sp_offset=spec[iw_min:iw_max],
            remove_mean=True,
            frac_cosbell=0.10,
            zero_padding=50,
            debugplot=0
        )

        offsets[i] = off
        fpeaks[i] = fp

    print("INFO: Cross-correlation completed")
    if debug_level >= 1:
        plot_offsets(offsets)
    
    return offsets, fpeaks


def crosscorrelate_spectra(cube_data, wavelength, config_parameters,
                           table_results_fitting, bin_map=None):
    """
    Cross-correlate spectra against a reference spectrum.
    The function uses "periodic_corr1d" from numina package.

    Works for:
    - pixel-by-pixel cubes (bin_map=None)
    - Voronoi-binned cubes (bin_map provided)

    Returns offset and correlation peak per spaxel.
    """

    # Get the reference spectrum from the YAML configuration file:
    x_ref, y_ref = np.array(config_parameters["ref_spec"])  # FITS
    x_ref -= 1
    y_ref -= 1

    if bin_map is None:
        spectrum_ref = cube_data[:, y_ref, x_ref]   # get the original spectrum at (xref, yref)
    else:
        ref_bin = bin_map[y_ref, x_ref]
        spectrum_ref = cube_data[:, ref_bin]        # get the spectrum of the Voronoi bin to which (xref, yref) belongs

    # The cross-correlation region is the whole continuum region without excluding anything
    lambda_min, lambda_max = config_parameters["reg_continuum"]
    nw = len(wavelength)

    iw_min = max(np.searchsorted(wavelength, lambda_min), 0)
    iw_max = min(np.searchsorted(wavelength, lambda_max), nw)

    spectrum_region_ref = spectrum_ref[iw_min:iw_max]

    # pixel-by-pixel cube (no voronoi binning done)
    if bin_map is None:
        nw, ny, nx = cube_data.shape
        #offsets = np.full((ny, nx), np.nan) # they will be saved as a 2D map
        #fpeaks  = np.full((ny, nx), np.nan) # they will be saved as a 2D map
        offsets_pixel = []
        fpeaks_croscorr = []

        for y in tqdm(range(ny), desc="Cross-correlating spectra pixel by pixel"):
            for x in range(nx):
                spec = cube_data[:, y, x]           # extract one spectrum
                if np.all(~np.isfinite(spec)):      # skip empty spectra
                    continue

                off, fp = periodic_corr1d(
                    sp_reference=spectrum_region_ref,
                    sp_offset=spec[iw_min:iw_max],
                    remove_mean=True,
                    frac_cosbell=0.10,
                    zero_padding=50,
                    debugplot=0   # 0: no plots; 12: one plot (crosscorrelation only); 22: more plots
                )

                #offsets[y, x] = off
                #fpeaks[y, x] = fp
                offsets_pixel.append(off)
                fpeaks_croscorr.append(fp)

        offsets_pixel_array = np.array(offsets_pixel)
        fpeak_croscorr_array = np.array(fpeaks_croscorr)

        return offsets_pixel_array, fpeak_croscorr_array  #offsets.ravel(), fpeaks.ravel() #.ravel() flattens an NDarray to 1D without changing the order
    
    # Voronoi-binned cube
    else:
        n_bins = cube_data.shape[1]
        offsets_bin = np.zeros(n_bins)
        fpeaks_croscorr_bin = np.zeros(n_bins)

        for i in range(n_bins):
            spec = cube_data[:, i]

            off, fp = periodic_corr1d(
                sp_reference=spectrum_region_ref,
                sp_offset=spec[iw_min:iw_max],
                remove_mean=True,
                frac_cosbell=0.10,
                zero_padding=50,
                debugplot=0
            )

            offsets_bin[i] = off                    # shape (nbins,) The offset for the bin ID "i"
            fpeaks_croscorr_bin[i] = fp             # shape (nbins,)

        # Propagate bin values back to spaxels (1D!)
        bin_id = table_results_fitting["bin_id"]    # shape (N_spaxels,) The bin ID of every spexel in the original cube
        offsets = offsets_bin[bin_id]
        fpeaks  = fpeaks_croscorr_bin[bin_id]
        print(len(offsets), len(table_results_fitting))

        return offsets, fpeaks



def crosscorrelate_spectra_datacube(cube_data, wavelength, config_parameters):

    nw, ny, nx = cube_data.shape

    x_ref, y_ref = np.array(config_parameters["ref_spec"]) # FITS format
    x_ref -= 1 # change to python indices
    y_ref -= 1
    spectrum_ref = cube_data[:, y_ref, x_ref] # correct to make it python-index

    crosscorr_region = np.array(config_parameters["reg_continuum"])
    print(f"The region in which spectra are being crosscorrelated is: {crosscorr_region} Å")

    lambda_min_crosscorr, lambda_max_crosscorr = crosscorr_region

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
    
    """
    Convert pixel offsets into velocity.

    Parameters
    ----------
    offsets : ndarray (N,)
    wavelength : ndarray (n_lambda,)
    """

    dlambda_dp = np.mean(np.diff(wavelength)) # Å/pixel o nm/pixel
    lambda_obs = np.array(line_restframe) * (1 + redshift)

    # Considering linear wavelength axis:
    delta_lambda = offsets_pixel_array * dlambda_dp
    c_kms = 299792.458
    velocity_array = c_kms * delta_lambda / lambda_obs
    print(f"INFO: The velocity values for every offset have been calculated: {velocity_array} km/s")
    print(f"INFO: Velocity computed for {len(velocity_array)} spectra")

    return velocity_array



