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
from ..model_fitting.spectrum_fitting import fit_gaussians_to_all_spectra_lmfit
from ..model_fitting.debug_plot_random_fits import debug_random_fits

def measure_spectra_properties(
    spectra,
    wavelength_range,
    config_parameters,
    analysis_table,
    redshift,
    line_restframe,
    debug_level
):
    """
    Apply physical measurements to spectra.
    Works for real AND simulated spectra.
    """
    print(f"INFO:Crosscorrelation to reference spectrum")
    offsets, fpeaks = crosscorrelate_spectra_unified(
        spectra,
        wavelength_range,
        config_parameters,
        analysis_table,
        debug_level
    )

    print(f"INFO: Calculating velocities from offsets")
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

    input("press enter to continue")
    if config_parameters["line_model_fit"]:
        print('INFO: Fitting line models with lmfit')
        if debug_level>=1:
            print("INFO: Generating debug plots for random spectra")
            debug_random_fits(
                spectra,
                wavelength_range,
                config_parameters,
                n_examples=10
            )

        analysis_table = fit_gaussians_to_all_spectra_lmfit(
            spectra,
            wavelength_range,
            analysis_table,
            config_parameters
        )
        analysis_table[:5].pprint()
        analysis_table[-5:].pprint()

    return analysis_table