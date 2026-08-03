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
        "mu",
        "sigma",
        "area_ha",
        "area_nii6548",
        "area_nii6583",
        "area_total",
        "chi2",
    ],

    "double_gaussian": [
        "amp1",
        "amp2",
        "mu1",
        "mu2",
        "sigma",
        "area1",
        "area2",
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
        "mu": "mu_gauss",
        "sigma": "sigma_gauss",
        "area_ha": "area_ha",
        "area_nii6548": "area_nii6548",
        "area_nii6583": "area_nii6583",
        "area_total": "area_total",
        "chi2": "chi2_gauss",
    },

    "double_gaussian": {
        "amp1": "amp1",
        "amp2": "amp2",
        "mu1": "mu1_gauss",
        "mu2": "mu2_gauss",
        "sigma": "sigma_gauss",
        "area1": "area1",
        "area2": "area2",
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
        "mu",
        "sigma",
        "area_ha",
        "area_nii6548",
        "area_nii6583",
        "area_total",
        "chi2",
    ],

    "double_gaussian": [
        "amp1",
        "amp2",
        "mu1",
        "mu2",
        "sigma",
        "area1",
        "area2",
        "area_total",
        "chi2",
    ]
}

