#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np 

def get_region_mask(
    wavelength,
    center=None,
    window=None,
    region=None):
    """
    Returns a boolean mask for a wavelength region.
    """
    if window is not None:
        if center is None:
            raise ValueError("Center wavelength required when using window.")
        left = center - window
        right = center + window
    else:
        left, right = region
    
    real_region_mask = (wavelength > left) & (wavelength < right)
    left_real=np.min(wavelength[real_region_mask])
    right_real=np.max(wavelength[real_region_mask])

    return real_region_mask, left_real, right_real


def apply_excluded_regions(mask, wavelength, excluded_regions):
    """
    Remove excluded wavelength regions from a mask.
    """
    if excluded_regions is None:
        return mask

    for ex_left, ex_right in excluded_regions:
        mask &= ~((wavelength > ex_left) & (wavelength < ex_right))

    return mask
