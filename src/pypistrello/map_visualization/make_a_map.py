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
import os
import teareduce as tea

from .build_2d_map import build_2d_map
from .calculate_contours import calculate_contours
from .calculate_contours import add_contours_to_plot


def make_a_map(table, wcs, config_parameters, working_dir, output_dir_path, map_choice):
    """
    Create a 2D map from tabular data and apply visualization settings.

    This function generates a specific map
    (e.g., flux, velocity, dispersion) using the input table and WCS.

    Parameters
    ----------
    table : astropy.table.Table
        Table containing the spatial coordinates and physical quantities
        derived from the spectral analysis.
    wcs : astropy.wcs.WCS
        World Coordinate System associated with the FITS data,
        used for proper spatial projection.
    config_parameters : dict
        Dictionary of parameters loaded from the YAML configuration file.
        These parameters control plotting style, contour settings,
        and map-specific options.
    working_dir : str
        Path to the working directory where temporary or auxiliary files
        (e.g., saved contours) will be written.
    output_dir_path : str
        Directory where the final map images will be saved.
    map_choice : str
        Identifier of the map being generated, used to select
        the appropriate configuration parameters."""
    
    # Prepare the parameters in the YAML in the correct format to be used
    if map_choice == "flux":
        yaml_key = "flux_map"
    elif map_choice == "vel":
        yaml_key = "velocity_map"
    else:
        raise ValueError("map_choice must be 'flux' or 'vel'")
    
    params = config_parameters[yaml_key]
    visualize = params.get("visualize", False)
    interpolate = params.get("interpolate", True)
    interp_method = params.get("interpolation_method", "nearest")
    zscale_factor = params.get("zscale_factor")
    if zscale_factor is None:
        zscale_factor = 0.05

    # Extract column information from the Table
    x = table["x"] # FITS format from the Table Attention!
    y = table["y"] # FITS format from the Table
    data = table[params["data_column"]]

    # Now we build a 2d map, the result will be different depending if interpolation is True:
    zi = build_2d_map(x, y, data, interpolate=interpolate, method=interp_method)

    # PLOT SET UP ----------------------------------
    
    fig = plt.figure(figsize=(7, 6))

    # Axis labels
    if params.get("wcs_activate", False) and wcs is not None:
        ax = plt.subplot(projection=wcs)
        ax.set_xlabel("RA")
        ax.set_ylabel("DEC")
    else:
        ax = plt.subplot()
        ax.set_xlabel("X [pix]")
        ax.set_ylabel("Y [pix]")

    # Min and max values for colorbar
    vmin = params.get("vmin")
    vmax = params.get("vmax")
    if vmin is None and vmax is None:
        finite_zi = zi[np.isfinite(zi)]
        if finite_zi.size == 0:
            raise ValueError("No finite data available for zscale")
        vmin, vmax = tea.zscale(image=finite_zi, factor=zscale_factor)
        print(f"Auto zscale applied: vmin={vmin:.3e}, vmax={vmax:.3e}")

    # Plotting
    im = ax.imshow(zi, origin="lower",
        cmap=params.get("cmap", "viridis"),
        vmin=vmin, vmax=vmax)

    # Colorbar settings
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(params.get("colorbar_label", ""))

    # Contours
    if params.get("calculate_contours", False): 
        print("Calculating contours")
        calculate_contours(params, table, working_dir, map_choice, ax)
        
    if params.get("add_flux_contours", False):
        contour_file_loaded = params.get("flux_contour_file")
        if contour_file_loaded is None:
            raise ValueError(
                "add_flux_contours=True but no flux_contour_file specified"
            )

        contour_file_loaded_path = working_dir / contour_file_loaded
        print(f"Adding contours from file {contour_file_loaded_path}")
        add_contours_to_plot(params, contour_file_loaded, ax)
        
    # Save the map
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

    print(f"Map saved in {output_path}")
    