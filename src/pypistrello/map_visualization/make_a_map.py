#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from cmap import Colormap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.pyplot as plt
import numpy as np
import os
import teareduce as tea

from .build_2d_map import build_2d_map
from .calculate_contours_flux import prepare_map_for_contours, compute_contour_levels, save_processed_map, draw_contours, load_processed_map


def make_a_map(table, wcs, config_parameters, working_dir,
               output_dir_path, map_choice):
    """
    Create a 2D map from tabular data and apply visualization settings.

    Works for:
    - spaxel-based tables (no Voronoi)
    - bin-based tables (Voronoi)

    Parameters
    ----------
    table : astropy.table.Table
        Table containing spatial coordinates and physical quantities.
        - If no Voronoi: one row per spaxel
        - If Voronoi: one row per bin
    wcs : astropy.wcs.WCS
    config_parameters : dict
    working_dir : str or Path
    output_dir_path : str or Path
    map_choice : str
    bin_map : ndarray (ny, nx), optional
        Required if table is Voronoi-binned
    """

    # Select YAML parameters
    if map_choice == "flux":
        yaml_key = "flux_map"
    elif map_choice == "vel":
        yaml_key = "velocity_map"
    elif map_choice == "snr":
        yaml_key = "snr_map"
    elif map_choice == "voronoi":
        yaml_key = "voronoi_map"
    elif map_choice == "sigma":
        yaml_key = "sigma_map"
    elif map_choice == "EW":
        yaml_key = "equivalent_width_map"
    else:
        raise ValueError("map_choice must be 'flux', 'vel', 'sigma', 'snr', 'voronoi' or 'EW'")

    params = config_parameters[yaml_key]

    visualize = params.get("visualize", False)
    interpolate = params.get("interpolate", True)
    interp_method = params.get("interpolation_method", "nearest")

    vcenter = params.get("vcenter")
    zscale_factor = params.get("zscale_factor", 0.05)

    cmap_name = params.get("cmap", "viridis")
    colorbar_loc = params.get("colorbar_loc", "right")
    print(f"[INFO] Colormap: {cmap_name}")
    cmap = Colormap(cmap_name).to_mpl()

    data_column = params["data_column"]

    # Build the map
    print("INFO: Building map from individual spaxels")
    x = table["x"]
    y = table["y"]
    data = table[data_column]

    zi = build_2d_map(
        x, y, data,
        interpolate=interpolate,
        method=interp_method
    )

    print(f"Map shape: {zi.shape}")

    # Plot setup
    fig = plt.figure(figsize=(7, 6))

    if params.get("wcs_activate", False) and wcs is not None:
        ax = plt.subplot(projection=wcs)
        ax.set_xlabel("RA")
        ax.set_ylabel("DEC")
    else:
        ax = plt.subplot()
        ax.set_xlabel("X [pix]")
        ax.set_ylabel("Y [pix]")

    ax.tick_params(
        direction="in",
        length=6,
        width=1.2,
    )

    # Color scaling
    vmin = params.get("vmin")
    vmax = params.get("vmax")

    if vmin is None and vmax is None:
        finite_zi = zi[np.isfinite(zi)]

        if finite_zi.size == 0:
            raise ValueError("No finite data available for zscale")

        vmin, vmax = tea.zscale(image=finite_zi, factor=zscale_factor)
        print(f"[INFO] Auto zscale: vmin={vmin:.3e}, vmax={vmax:.3e}")

    # Plot image
    if vcenter is not None:
        print(f"[INFO] Using TwoSlopeNorm centered at {vcenter}")
        norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
        im = ax.imshow(zi, origin="lower", cmap=cmap, norm=norm)
    else:
        im = ax.imshow(zi, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)

    # colorbar
    #cbar = plt.colorbar(im, ax=ax, pad=0.0, location=colorbar_loc)
    #cbar.set_label(params.get("colorbar_label", ""))

    if colorbar_loc == "right":

        cax = ax.inset_axes([1.00, 0.0, 0.05, 1.0], transform=ax.transAxes,)

        cbar = fig.colorbar(im, cax=cax, orientation="vertical",)
        cbar.set_label(params.get("colorbar_label", ""))

    elif colorbar_loc == "top":

        cax = ax.inset_axes([0.0, 1.00, 1.0, 0.05], transform=ax.transAxes,)

        cbar = fig.colorbar(im, cax=cax, orientation="horizontal",)
        cbar.set_label(params.get("colorbar_label", ""))

        cbar.ax.xaxis.set_ticks_position("top")
        cbar.ax.xaxis.set_label_position("top")

    else:
        raise ValueError("colorbar_loc must be 'right' or 'top'")

    
    if params.get("calculate_contours", False):

        contour_image, contour_mask = prepare_map_for_contours(zi, smooth_sigma=params.get("smooth_sigma", 1.0), fill_nan=True,)
        contour_levels = compute_contour_levels(contour_image, params,)

        save_processed_map(working_dir / f"{map_choice}_contours.npz", contour_image, contour_levels,)

        draw_contours(ax, contour_image, contour_levels,
            colors=params.get("contours_color", "black"),
            linewidths=params.get("contours_linewidth", 1.0),
            contours_outlined=params.get("contours_outlined", True)
        )

    if params.get("add_flux_contours", False):

        image, levels = load_processed_map(working_dir / params["flux_contour_file"])

        draw_contours(ax, image, levels,
            colors=params.get("contours_color", "black"),
            linewidths=params.get("contours_linewidth", 1.0),
            contours_outlined=params.get("contours_outlined", False)
        )

    # Save figure
    bg = params.get("background", None)

    if bg == "transparent":
        transparent = True
    else:
        transparent = False
        if bg is not None:
            fig.patch.set_facecolor(bg)
            ax.set_facecolor(bg)

    os.makedirs(output_dir_path, exist_ok=True)

    output_path = os.path.join(output_dir_path, f"{map_choice}_map.pdf")
    output_path_png = os.path.join(output_dir_path, f"{map_choice}_map.png")

    plt.savefig(output_path, dpi=200, bbox_inches="tight")

    if transparent:
        fig.savefig(output_path_png, dpi=300, transparent=True)

    if visualize:
        plt.show()

    plt.close(fig)

    print(f"[INFO] Map saved in {output_path}")