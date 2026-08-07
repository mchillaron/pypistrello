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
    ]
}

columns_to_copy = ["n_pix", "bin_center_x", "bin_center_y", 
                    "bin_area_trapz","bin_cont_noise","bin_snr_trapz","bin_cont_coeffs",
                    "velocity", "offsets", 
                    "amp_gauss", "mu_gauss", "sigma_gauss", "fwhm", "cont_gauss", "area_gauss", "chi2_gauss",
                    "amp_ha", "amp_nii6548", "amp_nii6583", "mu_nii6548", "mu_nii6583", "area_ha", "area_nii6548", "area_nii6583", "area_total",
                    "amp_line1", "amp_line2", "mu_line1", "mu_line2", "area_line1", "area_line2", "area_total", 
                    "sigmavel_gauss_AA", "sigmavel_gauss_AA_corr", "sigmavel_gauss_kms", "sigmavel_gauss_kms_corr",
                    "sigmavel_line1_AA", "sigmavel_line1_AA_corr", "sigmavel_line1_kms", "sigmavel_line1_kms_corr", 
                    "sigmavel_line2_AA", "sigmavel_line2_AA_corr", "sigmavel_line2_kms", "sigmavel_line2_kms_corr",
                    "sigmavel_ha_AA", "sigmavel_ha_AA_corr", "sigmavel_ha_kms", "sigmavel_ha_kms_corr",
                    "sigmavel_nii6548_AA", "sigmavel_nii6548_AA_corr", "sigmavel_nii6548_kms", "sigmavel_nii6548_kms_corr",
                    "sigmavel_nii6583_AA", "sigmavel_nii6583_AA_corr", "sigmavel_nii6583_kms", "sigmavel_nii6583_kms_corr"]


