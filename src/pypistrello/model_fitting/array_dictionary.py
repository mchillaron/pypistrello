#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

FIELDS = {
    "gaussian": [
        "amp",
        "mu",
        "sigma",
        "area",
        "chi2",
    ],

    "gaussian_area_fixed": [
        "amp",
        "mu",
        "sigma",
        "area",
        "chi2",
    ],

    "triplet_hanii": [
        "amp_ha",
        "amp_nii6548",
        "amp_nii6583",
        "mu_ha",
        "mu_nii6548",
        "mu_nii6583",
        "sigma",
        "area_ha",
        "area_nii6548",
        "area_nii6583",
        "area_total",
        "chi2",
    ],

    "double_gaussian": [
        "amp_line1",
        "amp_line2",
        "mu_line1",
        "mu_line2",
        "sigma",
        "area_line1",
        "area_line2",
        "area_total",
        "chi2",
    ],

    "narrow_broad": [
        "amp_narrow",
        "amp_broad",
        "mu_narrow",
        "mu_broad",
        "sigma",            # narrow sigma (kept as 'sigma' for fwhm/dispersion)
        "sigma_broad",
        "area_narrow",
        "area_broad",
        "area_total",
        "frac_broad",
        "frac_broad_err",
        "chi2",
        "bic",
    ],
}


COLUMN_NAMES = {

    "gaussian": {
        "amp": "amp_gauss",
        "mu": "mu_gauss",
        "sigma": "sigma_gauss",
        "area": "area_gauss",
        "chi2": "chi2_gauss",
    },

    "gaussian_area_fixed": {
        "amp": "amp_gauss",
        "mu": "mu_gauss",
        "sigma": "sigma_gauss",
        "area": "area_gauss",
        "chi2": "chi2_gauss",
    },

    "triplet_hanii": {
        "amp_ha": "amp_ha",
        "amp_nii6548": "amp_nii6548",
        "amp_nii6583": "amp_nii6583",
        "mu_ha": "mu_gauss",
        "mu_nii6548": "mu_nii6548",
        "mu_nii6583": "mu_nii6583",
        "sigma": "sigma_gauss",
        "area_ha": "area_ha",
        "area_nii6548": "area_nii6548",
        "area_nii6583": "area_nii6583",
        "area_total": "area_total",
        "chi2": "chi2_gauss",
    },

    "double_gaussian": {
        "amp_line1": "amp_line1",
        "amp_line2": "amp_line2",
        "mu_line1": "mu_line1",
        "mu_line2": "mu_line2",
        "sigma": "sigma_gauss",
        "area_line1": "area_line1",
        "area_line2": "area_line2",
        "area_total": "area_total",
        "chi2": "chi2_gauss",
    },

    "narrow_broad": {
        "amp_narrow": "amp_narrow",
        "amp_broad": "amp_broad",
        "mu_narrow": "mu_gauss",
        "mu_broad": "mu_broad",
        "sigma": "sigma_gauss",
        "sigma_broad": "sigma_broad",
        "area_narrow": "area_narrow",
        "area_broad": "area_broad",
        "area_total": "area_total",
        "frac_broad": "frac_broad",
        "frac_broad_err": "frac_broad_err",
        "chi2": "chi2_gauss",
        "bic": "bic_gauss",
    },
}

EMPTY_RESULTS = {

    "gaussian": [
        "amp",
        "mu",
        "sigma",
        "area",
        "chi2",
    ],

    "gaussian_area_fixed": [
        "amp",
        "mu",
        "sigma",
        "area",
        "chi2",
    ],

    "triplet_hanii": [
        "amp_ha",
        "amp_nii6548",
        "amp_nii6583",
        "mu_ha",
        "mu_nii6548",
        "mu_nii6583",
        "sigma",
        "area_ha",
        "area_nii6548",
        "area_nii6583",
        "area_total",
        "chi2",
    ],

    "double_gaussian": [
        "amp_line1",
        "amp_line2",
        "mu_line1",
        "mu_line2",
        "sigma",
        "area_line1",
        "area_line2",
        "area_total",
        "chi2",
    ],

    "narrow_broad": [
        "amp_narrow",
        "amp_broad",
        "mu_narrow",
        "mu_broad",
        "sigma",
        "sigma_broad",
        "area_narrow",
        "area_broad",
        "area_total",
        "frac_broad",
        "frac_broad_err",
        "chi2",
        "bic",
    ],
}

columns_to_copy = ["n_pix", "bin_center_x", "bin_center_y", 
                    "bin_area_trapz","bin_cont_noise","bin_snr_trapz","bin_cont_coeffs",
                    "velocity", "offsets", 
                    "amp_gauss", "mu_gauss", "sigma_gauss", "fwhm", "cont_gauss", "area_gauss", "chi2_gauss",
                    "amp_ha", "amp_nii6548", "amp_nii6583", "mu_nii6548", "mu_nii6583", "area_ha", "area_nii6548", "area_nii6583", "area_total", 
                    "amp_line1", "amp_line2", "mu_line1", "mu_line2", "area_line1", "area_line2", 
                    "amp_narrow", "amp_broad", "mu_broad", "sigma_broad", "fwhm_broad",
                    "area_narrow", "area_broad", "frac_broad", "frac_broad_err", "bic_gauss",
                    "sigmavel_gauss_AA", "sigmavel_gauss_AA_corr", "sigmavel_gauss_kms", "sigmavel_gauss_kms_corr",
                    "sigmavel_line1_AA", "sigmavel_line1_AA_corr", "sigmavel_line1_kms", "sigmavel_line1_kms_corr", 
                    "sigmavel_line2_AA", "sigmavel_line2_AA_corr", "sigmavel_line2_kms", "sigmavel_line2_kms_corr",
                    "sigmavel_ha_AA", "sigmavel_ha_AA_corr", "sigmavel_ha_kms", "sigmavel_ha_kms_corr",
                    "sigmavel_nii6548_AA", "sigmavel_nii6548_AA_corr", "sigmavel_nii6548_kms", "sigmavel_nii6548_kms_corr",
                    "sigmavel_nii6583_AA", "sigmavel_nii6583_AA_corr", "sigmavel_nii6583_kms", "sigmavel_nii6583_kms_corr"]


