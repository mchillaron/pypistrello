#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def signal_to_noise(wavelength, flux,
                    cont_mask, cont_fit_func,
                    line_mask, line_flux_trapz,
                    verbose=False):
    
    """
    Compute the signal-to-noise ratio (SNR) of an emission line
    from a one-dimensional spectrum.

    This function estimates the SNR of an emission line by comparing
    the integrated line flux to the propagated noise of the underlying
    continuum. The signal is defined as the continuum-subtracted line
    flux integrated over the wavelength range of the emission line.
    The noise is estimated from the root-mean-square (RMS) of the
    continuum residuals and propagated over the spectral width of
    the line following standard error propagation.

    The procedure consists of the following steps:

    1. The continuum is evaluated using a pre-fitted model
       (e.g. a low-order polynomial) over the wavelength regions
       free of line emission.

    2. The noise per spectral channel is estimated as the RMS of
       the continuum-subtracted flux in the continuum regions.

    3. The integrated noise over the emission line is computed by
       propagating the per-channel noise assuming independent
       Gaussian errors. For a line spanning N spectral channels
       of width Δλ, the total noise is:

           σ_line = σ_cont * sqrt(N) * Δλ

    4. The signal-to-noise ratio is then defined as:

           SNR = F_line / σ_line

       where F_line is the integrated line flux.

    Notes
    -----
    - The noise is assumed to be approximately constant over the
      wavelength range of the emission line.
    - The wavelength sampling is assumed to be uniform.
    - This method follows standard treatments of noise propagation
      in spectroscopy (e.g. Bevington & Robinson 2003) and is commonly
      used in integral-field spectroscopy analyses (e.g. Cappellari
      & Copin 2003; Husemann et al. 2013, A&A, 549, A87).

    Parameters
    ----------
    wavelength : array_like
        Wavelength array of the spectrum.
    flux : array_like
        Flux density array of the spectrum.
    cont_mask : array_like (bool)
        Boolean mask selecting continuum regions.
    cont_fit_func : callable
        Function representing the fitted continuum model.
    line_mask : array_like (bool)
        Boolean mask selecting the emission-line region.
    line_flux_trapz : float
        Integrated line flux computed using numerical integration.
    verbose : bool, optional
        If True, print diagnostic information.

    Returns
    -------
    snr : float
        Signal-to-noise ratio of the emission line.
    """

    signal = line_flux_trapz    # the signal of the line

    # To calculate the noise of the continuum flux, we first need to extract the continuum fit from "flux"
    lam_cont = wavelength[cont_mask]        # get the lambda range of the whole continuum region
    flux_cont = flux[cont_mask]             # get the signal of the whole continuum region
    cont_fit = cont_fit_func(lam_cont)      # evaluate the continuum fit in the continuum region
    flux_cont_sub = flux_cont - cont_fit    # Continuum-subtracted residuals

    # Noise calculation for 1 pixel
    noise_pix = np.nanstd(flux_cont_sub)
    print(f"Continuum RMS (with continuum subtraction): {noise_pix:.3e}")

    if verbose:
        noise_check = np.nanstd(flux_cont)
        print(f"Continuum RMS (with continuum subtraction): {noise_pix:.3e}")
        print(f"Continuum RMS (raw continuum): {noise_check:.3e}")

    # Since the line flux is integrated in N wavelength channels of width dλ=dl:
    N = line_mask.sum()
    if N <= 1:
        return np.nan
    dl = np.mean(np.diff(wavelength))
    print(f"The line flux has been integrated in N={N} wavelength channels of width dλ={dl:.3e}")

    noise = noise_pix*dl*np.sqrt(N)

    # SNR
    snr = signal/noise
    return snr