#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np
import matplotlib.pyplot as plt

def plot_trapz_spectrum(
    ax,
    wavelength,
    flux,
    lambda_line_sel,
    flux_line_sel,
    lambda_cont,
    flux_cont,
    cont_fit_func,
    line_left,
    line_right,
    line_obs,
    number_spectrum,
    excluded_region=None,
    zoom_limits=None,
    poly_order_cont=1,
    coords=None,
    y_pad=0.10,
):
    x, y = coords
    # --- Continuum fit
    cont_fit_cont = cont_fit_func(lambda_cont)
    cont_fit = cont_fit_func(lambda_line_sel)

    # --- Main curves
    ax.step(wavelength, flux, where="mid", color="k", lw=1, label="Spectrum")
    ax.plot(lambda_cont, cont_fit_cont, "darkorchid", ls="--", lw=1, label="Continuum fit")
    ax.scatter(lambda_cont, flux_cont, color="darkorchid", s=9, label="Continuum points")
    ax.scatter(lambda_line_sel, flux_line_sel, color="red", s=9, label="Line points")

    # --- Filled trapezoidal area
    lam_full = np.concatenate(([line_left], lambda_line_sel, [line_right]))
    cont_full = np.concatenate(([cont_fit[0]], cont_fit, [cont_fit[-1]]))
    flux_full = np.concatenate(([flux_line_sel[0]], flux_line_sel, [flux_line_sel[-1]]))

    ax.fill_between(
        lam_full, cont_full, flux_full,
        where=flux_full >= cont_full,
        color="darkseagreen", alpha=0.20, label="Positive area"
    )
    ax.fill_between(
        lam_full, cont_full, flux_full,
        where=flux_full < cont_full,
        color="orangered", alpha=0.25, label="Negative area"
    )

    # --- Guide lines
    ax.axvline(line_left, color="gray", ls=":", lw=0.8)
    ax.axvline(line_right, color="gray", ls=":", lw=0.8)
    ax.axvline(line_obs, color="gray", ls="--", lw=1, label="Line center")

    # --- Excluded region
    if excluded_region is not None:
        for ex_left, ex_right in excluded_region:
            ax.axvspan(ex_left, ex_right, color="gray", alpha=0.15)

    # --- Zoom
    if zoom_limits is not None:
        ax.set_xlim(*zoom_limits)

    # --- Autoscale Y
    mask = (wavelength > ax.get_xlim()[0]) & (wavelength < ax.get_xlim()[1])
    ymin, ymax = flux[mask].min(), flux[mask].max()
    ymargin = (ymax - ymin) * y_pad if ymax > ymin else 1
    ax.set_ylim(ymin - ymargin, ymax + ymargin)

    # --- Continuum equation
    p = cont_fit_func.coefficients
    if poly_order_cont == 0:
        eq = f"f(λ) = {p[0]:.2e}"
    elif poly_order_cont == 1:
        eq = f"f(λ) = {p[0]:.2e}·λ + {p[1]:.2e}"
    else:
        eq = f"Polynomial (order {poly_order_cont})"

    ax.text(
        0.98, 0.97, eq, transform=ax.transAxes,
        ha="right", va="top", fontsize=9,
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="gray")
    )

    # --- Labels
    ax.set_xlabel("Wavelength [Å]")
    ax.set_ylabel(r"Flux [erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$]")
    ax.set_title(f"{number_spectrum} - Spectrum at ({x}, {y})")
    ax.grid(alpha=0.25)
