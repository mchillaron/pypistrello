#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Create an interactive diagnostic plot for spectral analysis.

This is the initial step before performing the spectral line analysis.
The spectrum plotted is the integrated spectrum from a user-defined region
of the FITS table.
The plot is intended to help the user to decide on the parameters and it is interactive.
The upper panel shows the integrated spectrum with markers for the spectral lines of interest.
The lower panel shows a zoomed-in view of a zone chosen by the user by clicking on the upper panel.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
from matplotlib.gridspec import GridSpec

from .line_wavelength_dict import emission_lines
from .calculate_integrated_spectrum import calculate_integrated_spectrum


def plot_diagnostic_spectra(
    cube_data,
    wavelength_range,
    diagnostic_spectra,
    output_dir_path,
    redshift,
    line_restframe,
):
    """
    Create an interactive diagnostic plot for spectral analysis.

    The top panel shows the integrated spectrum over a user-defined spatial region.
    The bottom panel shows a zoomed-in view selected interactively from the top panel.

    Interactive features
    --------------------
    - Drag on top panel: zoom into wavelength range
    - r : reset zoom
    - a : define fit region (2 clicks)
    - c : define continuum region (2 clicks)
    - x : define excluded region (2 clicks)
    - u : undo last region
    - Esc : cancel current selection mode
    - Hover (bottom): wavelength / flux readout
    """

    y_pad = 0.1
    x1, x2, y1, y2 = diagnostic_spectra
    line_restframe = np.array(line_restframe)

    integrated_spectrum = calculate_integrated_spectrum(
        cube_data, wavelength_range, diagnostic_spectra
    )

    state = {
        "mode": None,           # a / c / x
        "clicks": [],
        "regions": {"a": [], "c": [], "x": []},
        "patches": [],
    }

    MODE_STYLE = {
        "a": dict(color="lightskyblue", label="Fit region"),
        "c": dict(color="wheat", label="Continuum region"),
        "x": dict(color="orange", label="Excluded region"),
    }

    fig = plt.figure(figsize=(9, 8))
    gs = GridSpec(
        3, 4,
        height_ratios=[1, 2, 0.05],
        width_ratios=[3, 3, 3, 1],
        hspace=0.4,
        wspace=0.3,
    )

    ax_top = fig.add_subplot(gs[0, :])
    ax_bot = fig.add_subplot(gs[1, :3])
    ax_leg = fig.add_subplot(gs[1, 3])
    ax_leg.axis("off")
    ax_top.grid(False)
    ax_bot.grid(False)

    # Legend / cheatsheet
    ax_leg.text(
        0.05, 0.95,
        "Controls\n"
        "──────────\n"
        "Drag (top): Zoom\n"
        "r : reset zoom\n"
        "a : fit region\n"
        "c : continuum region\n"
        "x : exclude region\n"
        "u : undo region\n"
        "Esc : cancel mode",
        va="top",
        fontsize=9,
    )

    #----------------------------------------------------------------------
    def plot_spectrum(ax, wave, flux, lw):
        ax.clear()
        ax.plot(wave, flux, color="black", lw=lw)
        ax.set_xlabel("Wavelength [Å]")
        ax.set_ylabel(r"$\mathrm{Flux\ [erg\ cm^{-2}\ s^{-1}\ \AA^{-1}]}$")
        ax.grid(False)

    def set_ylims(ax, flux):
        ymin, ymax = flux.min(), flux.max()
        pad = (ymax - ymin) * y_pad if ymax > ymin else 1.0
        ax.set_ylim(ymin - pad, ymax + pad)

    def plot_emission_lines(ax, xmin, xmax):
        ymax = ax.get_ylim()[1]
        for name, lam_rest in emission_lines.items():
            lam_obs = lam_rest * (1 + redshift)
            if xmin <= lam_obs <= xmax:
                ax.axvline(lam_obs, color="thistle", lw=1.4, linestyle='dotted', alpha=0.8)
                ax.text(
                    lam_obs + 0.1, ymax, name,
                    rotation=90, va="top", ha="center",
                    fontsize=8, color="darkslateblue",
                )
            
        # User selected wavelength to analyse
        lam_target = line_restframe * (1 + redshift)
        if xmin <= lam_target <= xmax:
            ax.axvline(
                lam_target,
                color="royalblue",
                lw=1.3,
                alpha=0.4,
                zorder=5,
            )


    # Initial top panel
    plot_spectrum(ax_top, wavelength_range, integrated_spectrum, lw=0.6)
    set_ylims(ax_top, integrated_spectrum)
    ax_top.set_title(
        f"Integrated spectrum from pixels ({x1},{y1}) - ({x2},{y2})"
    )
    plot_emission_lines(
        ax_top,
        wavelength_range.min(),
        wavelength_range.max(),
    )
    ax_top.plot(wavelength_range, integrated_spectrum, c='k', lw=0.6)
    xlim_full = (wavelength_range.min(), wavelength_range.max())
    ax_top.set_xlim(xlim_full)
    ax_top.set_autoscale_on(False)
    zoom_span = {"patch": None}

    # Zoom callback
    def onselect(xmin, xmax):
        if xmin == xmax:
            return

        if zoom_span["patch"] is not None:
            zoom_span["patch"].remove()

        zoom_span["patch"] = ax_top.axvspan(
            xmin, xmax, color="springgreen", alpha=0.05
        )

        # Bottom panel
        ax_bot.clear()
        state["patches"].clear()

        mask = (wavelength_range >= xmin) & (wavelength_range <= xmax)
        if not mask.any():
            fig.canvas.draw_idle()
            return

        zoom_wave = wavelength_range[mask]
        zoom_flux = integrated_spectrum[mask]

        plot_spectrum(ax_bot, zoom_wave, zoom_flux, lw=1.5)
        ax_bot.set_xlim(xmin, xmax)
        set_ylims(ax_bot, zoom_flux)
        plot_emission_lines(ax_bot, xmin, xmax)
        print(f"Bottom panel zoom: xmin = {xmin:.2f}, xmax = {xmax:.2f} Å")

        fig.canvas.draw_idle()

    span = SpanSelector(
        ax_top,
        onselect,
        "horizontal",
        useblit=False,
        interactive=True,
        props=dict(alpha=0.05, facecolor="springgreen"),
    )

    # Keyboard handler
    def on_key(event):
        if event.key == "r":
            ax_bot.clear()
            fig.canvas.draw_idle()

        elif event.key in ["a", "c", "x"]:
            state["mode"] = event.key
            state["clicks"].clear()

        elif event.key == "u" and state["patches"]:
            patch = state["patches"].pop()
            patch.remove()

        elif event.key == "escape":
            state["mode"] = None
            state["clicks"].clear()

        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("key_press_event", on_key)

    
    # Mouse click handler (bottom panel)
    def on_click(event):
        if event.inaxes != ax_bot or state["mode"] is None:
            return
        if event.xdata is None:
            return

        state["clicks"].append(event.xdata)
        style = MODE_STYLE[state["mode"]]

        if len(state["clicks"]) == 2:
            x1c, x2c = sorted(state["clicks"])
            print(f"Region '{style['label']}': {x1c:.2f},{x2c:.2f} Å")

            patch = ax_bot.axvspan(
                x1c, x2c, color=style["color"], alpha=0.3
            )
            state["regions"][state["mode"]].append((x1c, x2c))
            state["patches"].append(patch)
            state["clicks"].clear()
            state["mode"] = None

        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("button_press_event", on_click)

    plt.show()
    fig.savefig(output_dir_path / "diagnostic_plot.png", dpi=200)
    plt.close(fig)
