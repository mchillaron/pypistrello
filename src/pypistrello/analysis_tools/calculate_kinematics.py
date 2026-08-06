#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from astropy import units as u

import numpy as np

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



def calculate_dispersion(config_parameters, analysis_table):
    """
    Calculate the observed and instrumental-corrected velocity dispersion.

    For each fitted Gaussian component, the function computes:

        - Observed sigma in Angstrom
        - Observed sigma in km/s
        - Instrumental-corrected sigma in Angstrom
        - Instrumental-corrected sigma in km/s

    The results are appended to the analysis table.

    Parameters
    ----------
    config_parameters : dict
        YAML configuration dictionary.

    analysis_table : astropy.table.Table
        Output table from the Gaussian fitting.

    Returns
    -------
    analysis_table : astropy.table.Table
        Updated table including the new dispersion columns.
    """

    C = 299792.458  # km/s
    R = config_parameters["instr_resolution"]

    # Find Gaussian centres
    center_columns = sorted([c for c in analysis_table.colnames if c.startswith("mu")])

    if len(center_columns) == 0:
        raise ValueError("No Gaussian centre ('mu') columns found.")

    # Find sigma columns
    if "sigma" in analysis_table.colnames:
        sigma_columns = ["sigma"] * len(center_columns)         # Single sigma shared by every component

    else:
        sigma_columns = sorted([c for c in analysis_table.colnames if c.startswith("sigma")])

        if len(sigma_columns) != len(center_columns):
            raise ValueError(
                f"Found {len(center_columns)} centre columns "
                f"({center_columns}) but {len(sigma_columns)} sigma columns "
                f"({sigma_columns})."
            )

    # Instrumental sigma km/s units
    sigma_instr_velocity = C / (2.35482 * R)

    for i, (mu_col, sigma_col) in enumerate(
            zip(center_columns, sigma_columns), start=1):

        mu = analysis_table[mu_col].data
        sigma_lambda = analysis_table[sigma_col].data

        # Observed dispersion
        sigma_velocity = C * sigma_lambda / mu

        # Instrumental sigma in wavelength
        sigma_instr_lambda = mu / (2.35482 * R)

        # Corrected dispersions
        sigma_lambda_corr = np.sqrt(np.maximum(0.0, sigma_lambda**2 - sigma_instr_lambda**2))
        sigma_velocity_corr = np.sqrt(np.maximum(0.0, sigma_velocity**2 - sigma_instr_velocity**2))


        # Column suffix
        suffix = "" if len(center_columns) == 1 else str(i)

        analysis_table[f"sigmavel{suffix}_AA"] = (
            sigma_lambda * u.AA
        )

        analysis_table[f"sigmavel{suffix}_kms"] = (
            sigma_velocity * u.km / u.s
        )

        analysis_table[f"sigmavel{suffix}_AA_corr"] = (
            sigma_lambda_corr * u.AA
        )

        analysis_table[f"sigmavel{suffix}_kms_corr"] = (
            sigma_velocity_corr * u.km / u.s
        )

    return analysis_table