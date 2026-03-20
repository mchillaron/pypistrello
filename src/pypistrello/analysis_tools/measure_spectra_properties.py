#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import astropy.units as u
from ..area_fitting.crosscorrelation_spectra import crosscorrelate_spectra_unified
from ..area_fitting.crosscorrelation_spectra import convert_offset_velocity

def measure_spectra_properties(
    spectra,
    wavelength_range,
    config_parameters,
    analysis_table,
    redshift,
    line_restframe
):
    """
    Apply physical measurements to spectra.
    Works for real AND simulated spectra.
    """
    print(f"Crosscorrelation to reference spectrum")
    offsets, fpeaks = crosscorrelate_spectra_unified(
        spectra,
        wavelength_range,
        config_parameters,
        analysis_table
    )

    print(f"Calculating velocities from offsets")
    velocity = convert_offset_velocity(
        offsets,
        wavelength_range,
        redshift,
        line_restframe
    )

    analysis_table["offsets"] = offsets
    analysis_table["velocity"] = velocity * u.km / u.s
    analysis_table[:5].pprint()
    analysis_table[-5:].pprint()

    return analysis_table