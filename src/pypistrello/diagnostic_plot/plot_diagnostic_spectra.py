#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Function to create the diagnostic plot.
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

def calculate_integrated_spectrum(spectra_table, wavelength_range, diagnostic_spectra):
    x1, x2, y1, y2 = diagnostic_spectra

    integrated_spectrum = np.zeros_like(wavelength_range, dtype=float)
    n_spectra = 0
    for row in spectra_table:
        x = row["x"]
        y = row["y"]

        if x1 <= x <= x2 and y1 <= y <= y2:
            spectrum = row["spec"]
            if spectrum.shape != wavelength_range.shape:
                raise ValueError("Spectrum length does not match wavelength range.")
            
            integrated_spectrum += spectrum
            n_spectra += 1

    if n_spectra == 0:
        raise RuntimeError("No spectra found in selected region.")

    print(f"INFO: Integrated {n_spectra} spectra.")
    return integrated_spectrum



def plot_diagnostic_spectra(
        spectra_table,
        wavelength_range,
        diagnostic_spectra,
        output_dir_path,
        redshift,
        line_restframe,
    ):
        y_pad = 0.1
        x1, x2, y1, y2 = diagnostic_spectra

        integrated_spectrum = calculate_integrated_spectrum(
            spectra_table, wavelength_range, diagnostic_spectra)

        
        # Plot state
        state = {
            "plot_style": "line",   # "line" or "hist"
            "mode": None,           # f / c / x
            "clicks": [],
            "regions": {
                "f": [],
                "c": [],
                "x": [],
            },
        }

        MODE_STYLE = {
            "f": dict(color="cyan", label="Fit region"),
            "c": dict(color="purple", label="Continuum region"),
            "x": dict(color="orange", label="Excluded region"),
        }

        # Figure layout
        fig = plt.figure(figsize=(10, 12))
        gs = GridSpec(3, 4, height_ratios=[1, 2, 0.05], width_ratios=[3, 3, 3, 1])

        ax_top = fig.add_subplot(gs[0, :])
        ax_bot = fig.add_subplot(gs[1, :3])
        ax_leg = fig.add_subplot(gs[1, 3])
        ax_leg.axis("off")

        # Helper functions
        def plot_spectrum(ax, wave, flux, lw):
            ax.clear()
            if state["plot_style"] == "line":
                ax.plot(wave, flux, color="black", lw=lw)
            else:
                ax.step(wave, flux, where="mid", color="black", lw=lw)

            ax.set_xlabel("Wavelength [Å]")
            ax.set_ylabel(r"$\mathrm{Flux\ [erg\ cm^{-2}\ s^{-1}\ \AA^{-1}]}$")

        def set_ylims(ax, flux):
            ymin, ymax = flux.min(), flux.max()
            pad = (ymax - ymin) * y_pad if ymax > ymin else 1.0
            ax.set_ylim(ymin - pad, ymax + pad)

        def plot_emission_lines(ax, xmin, xmax):
            ymax = ax.get_ylim()[1]
            for name, lam_rest in emission_lines.items():
                lam_obs = lam_rest * (1 + redshift)
                if xmin <= lam_obs <= xmax:
                    ax.axvline(lam_obs, color="thistle", lw=1, alpha=0.8)
                    ax.text(
                        lam_obs + 0.1,
                        ymax,
                        name,
                        rotation=90,
                        va="top",
                        ha="center",
                        fontsize=8,
                        color="darkslateblue",
                    )

        # Initial top panel
        plot_spectrum(ax_top, wavelength_range, integrated_spectrum, lw=0.6)
        set_ylims(ax_top, integrated_spectrum)
        ax_top.set_title(
            f"Integrated spectrum in pixels ({x1},{y1}) - ({x2},{y2})"
        )
        plot_emission_lines(
            ax_top,
            wavelength_range.min(),
            wavelength_range.max(),
        )

        xlim_full = ax_top.get_xlim()

        
        # Zoom callback
        def onselect(xmin, xmax):
            print(f"Zoom selected: {xmin:.2f} – {xmax:.2f}")  
            
            mask = (wavelength_range >= xmin) & (wavelength_range <= xmax)
            if not mask.any():
                return

            zoom_wave = wavelength_range[mask]
            zoom_flux = integrated_spectrum[mask]

            plot_spectrum(ax_bot, zoom_wave, zoom_flux, lw=1.5)
            ax_bot.set_xlim(xmin, xmax)
            set_ylims(ax_bot, zoom_flux)
            plot_emission_lines(ax_bot, xmin, xmax)

            fig.canvas.draw_idle()

        span = SpanSelector(
            ax_top,
            onselect,
            "horizontal",
            useblit=True,
            props=dict(alpha=0.3, facecolor="thistle"),
        )

    
        # Reset zoom
        def reset_zoom():
            ax_bot.clear()
            ax_bot.set_xlabel("Wavelength [Å]")
            ax_bot.set_ylabel(r"$\mathrm{Flux\ [erg\ cm^{-2}\ s^{-1}\ \AA^{-1}]}$")
            fig.canvas.draw_idle()



        # Keyboard handler
        def on_key(event):
            if event.key == "k":
                state["plot_style"] = (
                    "hist" if state["plot_style"] == "line" else "line"
                )
                plot_spectrum(ax_top, wavelength_range, integrated_spectrum, lw=0.6)
                set_ylims(ax_top, integrated_spectrum)

                if ax_bot.has_data():
                    xmin, xmax = ax_bot.get_xlim()
                    onselect(xmin, xmax)

            elif event.key == "r":
                reset_zoom()

            elif event.key in ["f", "c", "x"]:
                state["mode"] = event.key
                state["clicks"].clear()
                print(f"Selection mode: {MODE_STYLE[event.key]['label']}")

            elif event.key == "escape":
                state["mode"] = None
                state["clicks"].clear()
                print("Selection cancelled")

            fig.canvas.draw_idle()

        fig.canvas.mpl_connect("key_press_event", on_key)

        # Mouse click handler
        def on_click(event):
            # Only bottom panel reacts to clicks
            if event.inaxes != ax_bot:
                return
            if state["mode"] is None:
                return

            x = event.xdata
            if x is None:
                return
            
            state["clicks"].append(x)
            style = MODE_STYLE[state["mode"]]

            ax_bot.axvline(x, color=style["color"], lw=1.5, alpha=0.9)

            if len(state["clicks"]) == 2:
                x1, x2 = sorted(state["clicks"])
                state["regions"][state["mode"]].append((x1, x2))
                print(f"{style['label']}: {x1:.2f} - {x2:.2f}")
                state["clicks"].clear()
                state["mode"] = None

            fig.canvas.draw_idle()

        fig.canvas.mpl_connect("button_press_event", on_click)

   
        # Hover readout
        annot = ax_bot.annotate(
            "",
            xy=(0, 0),
            xytext=(15, 15),
            textcoords="offset points",
            bbox=dict(boxstyle="round", fc="white"),
            arrowprops=dict(arrowstyle="->"),
        )
        annot.set_visible(False)

        def on_motion(event):
            if event.inaxes != ax_bot or event.xdata is None:
                annot.set_visible(False)
                fig.canvas.draw_idle()
                return

            annot.xy = (event.xdata, event.ydata)
            annot.set_text(
                f"λ = {event.xdata:.2f} Å\nF = {event.ydata:.2e}"
            )
            annot.set_visible(True)
            fig.canvas.draw_idle()

        fig.canvas.mpl_connect("motion_notify_event", on_motion)

        plt.show()

        # Save after interaction
        output_path = output_dir_path / "diagnostic_plot.png"
        fig.savefig(output_path, dpi=200)
        plt.close(fig)
