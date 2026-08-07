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
    Compute the observed and instrumental-corrected velocity dispersion
    for every Gaussian component found in the analysis table.

    New columns:
        sigmavel_<line>_AA
        sigmavel_<line>_AA_corr
        sigmavel_<line>_kms
        sigmavel_<line>_kms_corr
    """

    C = 299792.458          # km/s
    R = config_parameters["instr_resolution"]

    # Check Gaussian fit exists
    if "sigma_gauss" not in analysis_table.colnames:
        print("INFO: No Gaussian sigma found. Skipping velocity dispersion.")
        return analysis_table

    sigma_lambda = analysis_table["sigma_gauss"].data

    # Find all Gaussian areas (= physical lines)
    area_columns = sorted(
        c for c in analysis_table.colnames
        if c.startswith("area_")
        and c not in ("area_trapz", "area_total")
    )

    if len(area_columns) == 0:
        print("INFO: No Gaussian areas found.")
        return analysis_table

    sigma_instr_velocity = C / (2.35482 * R)

    # Loop over fitted lines
    for area_col in area_columns:

        line_name = area_col.replace("area_", "")

        if line_name == "ha":           # exception for Halpha+NII fitting
            mu_col = "mu_gauss"
        else:
            mu_col = f"mu_{line_name}"

        if mu_col not in analysis_table.colnames:
            print(f"WARNING: {mu_col} not found. Skipping.")
            continue

        mu = analysis_table[mu_col].data


        # Observed dispersion
        sigma_velocity = C * sigma_lambda / mu

        # Instrumental sigma
        sigma_instr_lambda = mu / (2.35482 * R)

        # Corrected dispersions
        sigma_lambda_corr = np.sqrt(np.maximum(0.0, sigma_lambda**2 - sigma_instr_lambda**2))
        sigma_velocity_corr = np.sqrt(np.maximum(0.0, sigma_velocity**2 - sigma_instr_velocity**2))

        # Column suffix
        analysis_table[f"sigmavel_{line_name}_AA"] = (sigma_lambda * u.AA)
        analysis_table[f"sigmavel_{line_name}_AA_corr"] = (sigma_lambda_corr * u.AA)
        analysis_table[f"sigmavel_{line_name}_kms"] = (sigma_velocity * u.km / u.s)
        analysis_table[f"sigmavel_{line_name}_kms_corr"] = (sigma_velocity_corr * u.km / u.s)

    return analysis_table


def calculate_equivalent_width(analysis_table):
    """
    Compute equivalent widths using both trapz measurements
    and Gaussian fits (if available).

    New columns:

        EW_trapz

    and, if Gaussian fits exist,

        EW_ha
        EW_nii6548
        EW_nii6583
        ...
    """

    # Continuum coefficients

    if "cont_coeffs" not in analysis_table.colnames:
        raise ValueError("Column 'cont_coeffs' not found.")

    coeffs_all = analysis_table["cont_coeffs"]

    # Trapz EW
    if "area_trapz" in analysis_table.colnames:

        ew = np.full(len(analysis_table), np.nan)

        # We evaluate the continuum at the Gaussian centre if available.
        # Otherwise we use the middle coefficient region (mean continuum).

        if "mu_gauss" in analysis_table.colnames:

            mu = analysis_table["mu_gauss"]

            for i, coeffs in enumerate(coeffs_all):
                cont = np.polyval(coeffs, mu[i])
                if cont > 0:
                    ew[i] = analysis_table["area_trapz"][i] / cont

        analysis_table["EW_trapz"] = ew * u.AA

    # Gaussian EW

    area_columns = sorted(
        c for c in analysis_table.colnames
        if c.startswith("area_")
        and c not in ("area_trapz", "area_total")
    )

    for area_col in area_columns:

        line_name = area_col.replace("area_", "")

        if line_name == "ha":
            mu_col = "mu_gauss"
        else:
            mu_col = f"mu_{line_name}"

        if mu_col not in analysis_table.colnames:
            continue

        mu = analysis_table[mu_col]
        ew = np.full(len(analysis_table), np.nan)

        for i, coeffs in enumerate(coeffs_all):

            continuum = np.polyval(coeffs, mu[i])
            if continuum > 0:
                ew[i] = analysis_table[area_col][i] / continuum

        analysis_table[f"EW_{line_name}"] = ew * u.AA

    return analysis_table


