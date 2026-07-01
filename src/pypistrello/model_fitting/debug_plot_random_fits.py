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
from tqdm import tqdm

from .fit_continuum_model import fit_continuum
from .line_models import fit_model_lmfit
from .spectrum_fitting import get_model_and_initial_params


def debug_random_fits(
    spectra,
    wavelength,
    analysis_table,
    config,
    offsets,
    n_examples=10
):
    """
    Plot random spectra with continuum + model fit for debugging.

    Supports:
      - gaussian
      - gaussian_area_fixed
      - triplet_hanii
    """

    print(f"INFO: Generating {n_examples} debug plots")

    idx_random = np.random.choice(spectra.shape[1], size=min(n_examples, spectra.shape[1]), replace=False)

    zoom = config.get("zoom_plot", [wavelength.min(), wavelength.max()])
    ypad = config.get("y_padding", 0.1)

    model_name = config.get(
        "model_to_fit",
        config.get("model", "gaussian")
    ).lower()

    for i in idx_random:

        spec = spectra[:, i]

        continuum_model = fit_continuum(wavelength, spec, config)
        if continuum_model is None:
            continue

        continuum = continuum_model(wavelength)
        spec_sub = spec - continuum

        # Fit line
        model_name = config.get("model", "gaussian").lower()

        if model_name == "triplet_hanii":
            lmin, lmax = config.get("reg_fitting_gaussians", config["reg_fitting"])
        elif model_name == "double_gaussian":
            lmin, lmax = config.get("reg_fitting_gaussians", config["reg_fitting"])
        else:
            lmin, lmax = config["reg_fitting"]

        mask = (wavelength >= lmin) & (wavelength <= lmax)
        x = wavelength[mask]
        y = spec_sub[mask]

        valid = np.isfinite(x) & np.isfinite(y)
        x = x[valid]
        y = y[valid]

        model_func, p0 = get_model_and_initial_params(config, x, y, wavelength, i, analysis_table)
        result = fit_model_lmfit(x, y, model_func, p0, config)

        if result is None:
            continue

        # High-resolution grid for smooth plotting
        if model_name == "triplet_hanii":

            center = result.params["center"].value
            sigma = result.params["sigma"].value

            HA_REST = 6562.80
            NII6548_REST = 6548.05
            NII6583_REST = 6583.45

            center_6548 = (center- (HA_REST - NII6548_REST))

            center_6583 = (center+ (NII6583_REST - HA_REST))

            x_plot_min = center_6548 - 5 * sigma
            x_plot_max = center_6583 + 5 * sigma

            x_fine = np.linspace(x_plot_min, x_plot_max, 2000)

        elif model_name=="double_gaussian":

            lmin_plot, lmax_plot = config.get(
                "reg_fitting_gaussians",
                config["reg_fitting"]
            )

            x_fine = np.linspace(lmin_plot, lmax_plot, 2000)

        else:
            x_fine = np.linspace(x.min(), x.max(), 1000)

        y_fine = result.eval(x=x_fine)               # Evaluate best-fit model on fine grid

        # Residuals in the fitting region:
        model_full = np.zeros_like(wavelength)
        model_mask = result.eval(x=wavelength[mask])
        model_full[mask] = model_mask
        spec_clean = spec_sub - model_full

        residuals = y - result.eval(x=x)
        rms = np.sqrt(np.mean(residuals**2))

        # Plot
        fig, (ax1, ax2) = plt.subplots(
            2, 1,
            figsize=(8, 6),
            sharex=True,
            gridspec_kw={"height_ratios": [3, 1]}
        )

        # Top panel: Spectrum + fit
        ax1.plot(wavelength, spec, color="black", alpha=0.5, label="Spectrum")
        ax1.plot(wavelength, continuum, "--", label="Continuum")

        ax1.plot(
            x_fine,
            y_fine + continuum_model(x_fine),
            linewidth=2,
            alpha=0.8,
            label="Model fit"
        )

        if model_name == "triplet_hanii":

            area_ha = result.params["area"].value
            amp_ha = area_ha / (sigma * np.sqrt(2*np.pi))
            amp_nii6583 = result.params["amp_nii6583"].value

            g_ha = amp_ha * np.exp(-0.5 * ((x_fine - center)/ sigma)**2)
            g_6583 = amp_nii6583 * np.exp(-0.5 * ((x_fine - center_6583)/ sigma)**2)
            g_6548 = (amp_nii6583/3.0) * np.exp(-0.5 * ((x_fine - center_6548)/ sigma)**2)

            cont_fine = continuum_model(x_fine)

            ax1.plot(
                x_fine,
                g_ha + cont_fine,
                "--",
                alpha=0.8,
                label="Hα"
            )

            ax1.plot(
                x_fine,
                g_6548 + cont_fine,
                "--",
                alpha=0.8,
                label="[NII]6548"
            )

            ax1.plot(
                x_fine,
                g_6583 + cont_fine,
                "--",
                alpha=0.8,
                label="[NII]6583"
            )
        
        elif model_name=="double_gaussian":

            sigma=result.params["sigma"].value
            center=result.params["center"].value
            delta=result.params["delta_lambda"].value

            center2=center+delta
            print(delta)
            print(center2)

            area1=result.params["area"].value
            print(area1)
            amp1=area1/(sigma*np.sqrt(2*np.pi))
            print(amp1)
            amp2=result.params["amp2"].value
            print(amp2)

            g1=amp1*np.exp(-0.5*((x_fine-center)/sigma)**2)
            g2=amp2*np.exp(-0.5*((x_fine-center2)/sigma)**2)

            cont=continuum_model(x_fine)

            ax1.plot(
                x_fine,
                g1+cont,
                "--",
                lw=1.8,
                label="Line 1"
            )

            ax1.plot(
                x_fine,
                g2+cont,
                "--",
                lw=1.8,
                label="Line 2"
            )
        
        ax1.set_xlim(zoom)
        
        ymin = np.nanmin(spec[(wavelength >= zoom[0]) & (wavelength <= zoom[1])])
        ymax = np.nanmax(spec[(wavelength >= zoom[0]) & (wavelength <= zoom[1])])
        dy = ymax - ymin

        ax1.set_ylim(ymin - ypad * dy, ymax + ypad * dy)

        ax1.tick_params(direction='in', which='both', top=True, right=True)
        ax1.minorticks_on()

        ax1.set_title(f"Spectrum {i}")
        ax1.legend()

        # Bottom panel: Residuals 
        ax2.axhline(0, linestyle="--")  # base line
        ax2.plot(
            wavelength,
            spec_clean,
            alpha=0.4,
            label="Spectrum - model"
        )

        # Residuals in fitting region
        ax2.plot(
            x,
            residuals,
            linewidth=1.7,
            color="hotpink",
            label="Fit residuals"
        )

        ax2.axvspan(lmin, lmax, alpha=0.1)

        ax2.set_xlim(zoom)
        ax2.set_ylabel("Residuals")
        ax2.set_xlabel("Wavelength")
        ax2.set_title(f"RMS = {rms:.3e}")

        # Automatic y-limits for residuals
        rmin = np.nanmin(residuals)
        rmax = np.nanmax(residuals)
        rdy = rmax - rmin if rmax > rmin else 1.0

        ax2.set_ylim(rmin - 0.1 * rdy, rmax + 0.1 * rdy)

        ax2.tick_params(direction='in', which='both', top=True, right=True)
        ax2.minorticks_on()

        plt.tight_layout()
        plt.show()







from matplotlib.backends.backend_pdf import PdfPages

def save_all_fits_to_pdf(
    spectra,
    wavelength,
    analysis_table,
    config,
    offsets,
    output_pdf="debug_spectra.pdf",
    ncols=2,
    nrows=2   
):
    """
    Save all spectra fits into a multi-page PDF with progress bar.
    Each plot contains spectrum + residuals.
    """

    n_spectra = spectra.shape[1]
    plots_per_page = ncols * nrows

    zoom = config.get("zoom_plot", [wavelength.min(), wavelength.max()])
    ypad = config.get("y_padding", 0.1)

    with PdfPages(output_pdf) as pdf:

        for start in tqdm(range(0, n_spectra, plots_per_page), desc="Generating PDF"):

            fig = plt.figure(figsize=(12, 4 * nrows))
            for j, i in enumerate(range(start, min(start + plots_per_page, n_spectra))):
                # GridSpec para subplot doble
                gs = fig.add_gridspec(
                    nrows * 3, ncols,
                    height_ratios=([4, 1, 0.8] * nrows),
                    hspace=0.30 #vertical space between ax1 and ax2 plots
                )

                row = j // ncols
                col = j % ncols

                base = row * 3

                #ax1 = fig.add_subplot(gs[2 * row, col])
                #ax2 = fig.add_subplot(gs[2 * row + 1, col], sharex=ax1)
                ax1 = fig.add_subplot(gs[base, col])
                ax2 = fig.add_subplot(gs[base + 1, col], sharex=ax1)

                spec = spectra[:, i]

                continuum_model = fit_continuum(wavelength, spec, config)
                if continuum_model is None:
                    continue

                continuum = continuum_model(wavelength)
                spec_sub = spec - continuum

                lmin, lmax = config["reg_fitting"]
                mask = (wavelength >= lmin) & (wavelength <= lmax)

                x = wavelength[mask]
                y = spec_sub[mask]

                model_func, p0 = get_model_and_initial_params(
                    config, x, y, wavelength, i, analysis_table
                )

                result = fit_model_lmfit(x, y, model_func, p0, config)

                if result is None:
                    continue
                
                # High-resolution grid for smooth plotting
                x_fine = np.linspace(x.min(), x.max(), 1000)
                y_fine = result.eval(x=x_fine)               # Evaluate best-fit model on fine grid

                # Complete model
                model_full = np.zeros_like(wavelength)
                model_full[mask] = result.eval(x=wavelength[mask])

                spec_clean = spec_sub - model_full
                residuals = y - result.eval(x=x)
                rms = np.sqrt(np.mean(residuals**2))

                # Main plot
                # Top panel: Spectrum + fit
                ax1.plot(wavelength, spec, color="black", alpha=0.5, label="Spectrum")
                ax1.plot(wavelength, continuum, "--", label="Continuum")

                ax1.plot(
                    x_fine,
                    y_fine + continuum_model(x_fine),
                    alpha=0.7,
                    label="Gaussian fit"
                )

                ax1.set_xlim(zoom)

                ymin = np.nanmin(spec[(wavelength >= zoom[0]) & (wavelength <= zoom[1])])
                ymax = np.nanmax(spec[(wavelength >= zoom[0]) & (wavelength <= zoom[1])])
                dy = ymax - ymin

                ax1.set_ylim(ymin - ypad * dy, ymax + ypad * dy)

                ax1.tick_params(direction='in', which='both', top=True, right=True)
                ax1.tick_params(labelbottom=False) # in case we want to hide x-axis labels in the top panel
                ax1.minorticks_on()

                ax1.set_title(f"Spectrum from Bin {i}")
                ax1.legend()

                # Bottom panel: Residuals 
                ax2.axhline(0, linestyle="--")  # base line
                ax2.plot(
                    wavelength,
                    spec_clean,
                    alpha=0.4,
                    label="Spectrum - model (full)"
                )

                # Residuals in fitting region
                ax2.plot(
                    x,
                    residuals,
                    linewidth=2,
                    color="hotpink",
                    label="Fit residuals"
                )

                ax2.axvspan(lmin, lmax, alpha=0.1)

                ax2.set_xlim(zoom)
                ax2.set_ylabel("Residuals")
                ax2.set_xlabel("Wavelength")
                ax2.set_title(f"RMS = {rms:.3e}")

                # Automatic y-limits for residuals
                rmin = np.nanmin(residuals)
                rmax = np.nanmax(residuals)
                rdy = rmax - rmin if rmax > rmin else 1.0

                ax2.set_ylim(rmin - 0.1 * rdy, rmax + 0.1 * rdy)

                ax2.tick_params(direction='in', which='both', top=True, right=True)
                ax2.minorticks_on()

            #plt.tight_layout()
            plt.tight_layout(pad=0.6)
            pdf.savefig(fig)
            plt.close(fig)

    print(f"PDF saved to {output_pdf}")



def save_log_spectra_pdf(
    spectra,
    wavelength,
    analysis_table,
    config,
    offsets,
    output_pdf="log_spectra.pdf",
    ncols=3,
    nrows=4,
    pedestal=0.1
):
    """
    Save spectra + Gaussian fits in log scale with pedestal.
    
    log(data + pedestal)
    """

    n_spectra = spectra.shape[1]
    plots_per_page = ncols * nrows

    zoom = config.get("zoom_plot", [wavelength.min(), wavelength.max()])
    ypad = config.get("y_padding", 0.1)

    with PdfPages(output_pdf) as pdf:

        for start in tqdm(range(0, n_spectra, plots_per_page), desc="Generating LOG PDF"):

            fig, axes = plt.subplots(
                nrows, ncols,
                figsize=(12, 10)
            )

            axes = axes.flatten()

            for j, i in enumerate(range(start, min(start + plots_per_page, n_spectra))):

                ax = axes[j]
                spec = spectra[:, i]

                continuum_model = fit_continuum(wavelength, spec, config)
                if continuum_model is None:
                    continue

                continuum = continuum_model(wavelength)
                spec_sub = spec - continuum

                # Región de ajuste
                lmin, lmax = config["reg_fitting"]
                mask = (wavelength >= lmin) & (wavelength <= lmax)

                x = wavelength[mask]
                y = spec_sub[mask]

                model_func, p0 = get_model_and_initial_params(
                    config, x, y, wavelength, i, analysis_table
                )
                result = fit_model_lmfit(x, y, model_func, p0, config)
                if result is None:
                    continue

                x_fine = np.linspace(x.min(), x.max(), 500)
                model_fine = result.eval(x=x_fine)

                # Complete model
                model_full = np.zeros_like(wavelength)
                model_full[mask] = result.eval(x=wavelength[mask])

                model_total = model_full + continuum

                # LOG + pedestal
                log_spec = np.log(spec + pedestal)
                log_cont = np.log(continuum + pedestal)
                log_model_fine = np.log(model_fine + continuum_model(x_fine) + pedestal)

                # PLOT
                ax.plot(wavelength, log_spec, alpha=0.4)
                ax.plot(wavelength, log_cont, "--")
                ax.plot(x_fine, log_model_fine)

                ax.axvspan(lmin, lmax, alpha=0.1)

                # Zoom
                ax.set_xlim(zoom)

                # Y-limits en log
                zoom_mask = (wavelength >= zoom[0]) & (wavelength <= zoom[1])
                yvals = log_spec[zoom_mask]

                ymin = np.nanmin(yvals)
                ymax = np.nanmax(yvals)
                dy = ymax - ymin if ymax > ymin else 1.0

                ax.set_ylim(ymin - ypad * dy, ymax + ypad * dy)

                ax.set_title(f"Spec {i}", fontsize=8)
                ax.tick_params(labelsize=6)

            # Ocultar ejes vacíos
            for k in range(j + 1, plots_per_page):
                axes[k].axis("off")

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    print(f"Log PDF saved to {output_pdf}")