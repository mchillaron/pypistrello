#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from .plotting import plot_trapz_spectrum

def save_trapz_plots_to_pdf(total_spectra, plot_inputs, pdf_path,
                            nrows=4, ncols=2):
    
    plots_per_page = nrows * ncols

    with PdfPages(pdf_path) as pdf:
        for page_start in range(0, len(plot_inputs), plots_per_page):

            fig, axes = plt.subplots(
                nrows, ncols,
                figsize=(8.5, 11),
                constrained_layout=True
            )

            axes = axes.flatten()

            #for ax, plot_kwargs in zip(
            #    axes,
            #    plot_inputs[page_start:page_start + plots_per_page]
            #):
            #    plot_trapz_spectrum(ax=ax, number_spectrum=j, **plot_kwargs)

            for i, (ax, plot_kwargs) in enumerate(
                zip(
                    axes,
                    plot_inputs[page_start:page_start + plots_per_page]
                )
            ):
                j = page_start + i + 1  # ⭐ GLOBAL spectrum index

                plot_trapz_spectrum(
                    ax=ax,
                    number_spectrum=j,
                    **plot_kwargs,
                )

            # Turn off unused axes
            for ax in axes[len(plot_inputs[page_start:page_start + plots_per_page]):]:
                ax.axis("off")

            pdf.savefig(fig)
            plt.close(fig)
