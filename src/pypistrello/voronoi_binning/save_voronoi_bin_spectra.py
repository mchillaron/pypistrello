#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

def save_voronoi_bin_spectra_pdf(
    spectra,
    wavelength,
    analysis_table,
    config_parameters,
    output_dir,
    simulated_spectra=None,
    grid_size=(5, 5),
    sort_by_snr=True
):
    """
    Save Voronoi bin spectra in a multi-page PDF (grid per page).

    Parameters
    ----------
    spectra : ndarray (n_lambda, n_bins)
    wavelength : ndarray (n_lambda,)
    analysis_table : Table (n_bins)
        Must contain 'bin_id', optionally 'bin_snr'
    simulated_spectra : ndarray (n_lambda, n_bins, Nsim), optional
    grid_size : tuple (nrows, ncols)
    sort_by_snr : bool
    """

    # Output path
    pdf_path = output_dir / "voronoi_bin_spectra.pdf"

    print(f"INFO: Saving PDF with summed spectra: {pdf_path}")

    # Plot settings
    wmin, wmax = config_parameters["zoom_plot"]
    iw_min = np.searchsorted(wavelength, wmin)
    iw_max = np.searchsorted(wavelength, wmax)
    wl = wavelength[iw_min:iw_max]

    # Sorting
    n_bins = spectra.shape[1]
    indices = np.arange(n_bins)

    if sort_by_snr and "bin_snr" in analysis_table.colnames:
        print("INFO: Sorting bins by SNR (descending)")
        indices = np.argsort(analysis_table["bin_snr"])[::-1]

    # PDF grid configuration
    nrows, ncols = grid_size
    per_page = nrows * ncols

    print(f"INFO: Grid: {nrows}x{ncols} → {per_page} bins/page")

    # PDF 
    with PdfPages(pdf_path) as pdf:

        for start in range(0, n_bins, per_page):

            fig, axes = plt.subplots(
                nrows, ncols,
                figsize=(ncols * 3, nrows * 2.5)
            )

            axes = axes.flatten()

            subset = indices[start:start + per_page]

            for ax, idx in zip(axes, subset):

                spec = spectra[iw_min:iw_max, idx]

                if np.all(~np.isfinite(spec)):
                    ax.axis("off")
                    continue

                # Histogram-style spectrum
                ax.step(wl, spec, where="mid", linewidth=1)

                # Uncertainty (simulations)
                if simulated_spectra is not None:
                    sim = simulated_spectra[iw_min:iw_max, idx, :]
                    mean = np.nanmean(sim, axis=1)
                    std = np.nanstd(sim, axis=1)

                    ax.fill_between(
                        wl,
                        mean - std,
                        mean + std,
                        alpha=0.3,
                        step="mid"
                    )

                # titles
                bin_id = analysis_table["bin_id"][idx]

                if "bin_snr" in analysis_table.colnames:
                    snr = analysis_table["bin_snr"][idx]
                    title = f"{bin_id} | SNR={snr:.1f}"
                else:
                    title = f"{bin_id}"

                ax.set_title(title, fontsize=7)

                # Regions
                if "reg_fitting" in config_parameters:
                    l1, l2 = config_parameters["reg_fitting"]
                    ax.axvspan(l1, l2, alpha=0.1)

                if "reg_continuum" in config_parameters:
                    l1, l2 = config_parameters["reg_continuum"]
                    ax.axvspan(l1, l2, alpha=0.05)

                # Minimal ticks
                ax.tick_params(labelsize=6)

            # Turn off empty panels
            for ax in axes[len(subset):]:
                ax.axis("off")

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            print(f"INFO: Saved page {start // per_page}")

    print("INFO: PDF completed")