#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import matplotlib.pyplot as plt
import numpy as np
from .fit_continuum_model import fit_continuum
from .line_models import fit_model_lmfit
from .spectrum_fitting import get_model_and_initial_params


def debug_random_fits(
    spectra,
    wavelength,
    config,
    n_examples=10
):
    """
    Plot random spectra with continuum + Gaussian fit for debugging.
    """

    print(f"INFO: Generating {n_examples} debug plots")

    idx_random = np.random.choice(spectra.shape[1], size=min(n_examples, spectra.shape[1]), replace=False)

    zoom = config.get("zoom_plot", [wavelength.min(), wavelength.max()])
    ypad = config.get("y_padding", 0.1)

    for i in idx_random:

        spec = spectra[:, i]

        continuum_model = fit_continuum(wavelength, spec, config)
        if continuum_model is None:
            continue

        continuum = continuum_model(wavelength)
        spec_sub = spec - continuum

        # Fit line
        lmin, lmax = config["reg_fitting"]
        mask = (wavelength >= lmin) & (wavelength <= lmax)

        x = wavelength[mask]
        y = spec_sub[mask]

        model_func, p0 = get_model_and_initial_params(config, x, y)
        result = fit_model_lmfit(x, y, model_func, p0)

        # High-resolution grid for smooth plotting
        x_fine = np.linspace(x.min(), x.max(), 1000)
        y_fine = result.eval(x=x_fine)               # Evaluate best-fit model on fine grid

        # Plot
        plt.figure(figsize=(8, 5))

        # Full spectrum
        plt.plot(wavelength, spec, color="black", alpha=0.5, label="Spectrum")

        # Continuum
        plt.plot(wavelength, continuum, "--", label="Continuum")

        # Fit
        if result is not None:
            #plt.plot(x, result.best_fit + continuum[mask], label="Gaussian fit")
            plt.plot(x_fine, y_fine + continuum_model(x_fine), label="Gaussian fit")

        # Zoom
        plt.xlim(zoom)

        ymin = np.nanmin(spec[(wavelength >= zoom[0]) & (wavelength <= zoom[1])])
        ymax = np.nanmax(spec[(wavelength >= zoom[0]) & (wavelength <= zoom[1])])
        dy = ymax - ymin

        plt.ylim(ymin - ypad * dy, ymax + ypad * dy)

        plt.title(f"Spectrum {i}")
        plt.legend()
        plt.show()