#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from pathlib import Path

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import teareduce as tea

from .file_loading.load_yaml_file import load_yaml_file
from .file_loading.load_fits_table import load_fits_table
from .line_fitting.map_plots import build_2d_map
from .line_fitting.map_plots import save_contours
from .line_fitting.map_plots import load_contours

RED     = "\033[91m"
GREEN   = "\033[92m"
YELLOW  = "\033[93m"
BLUE    = "\033[94m"
MAGENTA = "\033[95m"
CYAN    = "\033[96m"
BOLD = "\033[1m"
RESET   = "\033[0m"

def map_plotting(working_dir, fits_path, config_path, output_dir_path, map_choice):
    """Plotting maps from spectral line analysis results
    Parameters
    ----------
    """
    
    # Extract parameters from YAML configuration file:
    print(f"{BLUE}{BOLD} Reading parameters from YAML file{RESET}")
    config_parameters = load_yaml_file(config_path)
    print("YAML file read successfully")

    # Prepare the parameters in the YAML in the correct format to be used
    if map_choice == "flux":
        yaml_key = "flux_map"
    elif map_choice == "vel":
        yaml_key = "velocity_map"
    else:
        raise ValueError("map_choice must be 'flux' or 'vel'")
    
    params = config_parameters[yaml_key]
    interpolate = params.get("interpolate", True)
    interp_method = params.get("interpolation_method", "nearest")
    
    # read the table
    table, wcs = load_fits_table(fits_path)

    x = table["x"]
    y = table["y"]
    data = table[params["data_column"]]
    
    # Build a 2D map
    #xi, yi, zi = build_2d_map(x, y, data)
    zi = build_2d_map(
        x,
        y,
        data,
        interpolate=interpolate,
        method=interp_method
    )

    #Plot setup
    fig = plt.figure(figsize=(7, 6))

    if params.get("wcs_activate", False) and wcs is not None:
        ax = plt.subplot(projection=wcs)
        ax.set_xlabel("RA")
        ax.set_ylabel("DEC")
    else:
        ax = plt.subplot()
        ax.set_xlabel("X [pix]")
        ax.set_ylabel("Y [pix]")

    vmin = params.get("vmin")
    vmax = params.get("vmax")

    if vmin is None and vmax is None:
        finite_zi = zi[np.isfinite(zi)]
        if finite_zi.size == 0:
            raise ValueError("No finite data available for zscale")
        vmin, vmax = tea.zscale(image=finite_zi, factor=0.05)
        #vmin, vmax = tea.zscale(image=zi, factor=0.05)
        print(f"Auto zscale applied: vmin={vmin:.3e}, vmax={vmax:.3e}")

    im = ax.imshow(
        zi,
        origin="lower",
        cmap=params.get("cmap", "viridis"),
        vmin=vmin,
        vmax=vmax
    )

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(params.get("colorbar_label", ""))
    
    # contours
    if params.get("calculate_contours", False):
        zi_contours = build_2d_map(
            x, y, data,
            interpolate=True,
            method="linear"
        )

        levels = params.get("contour_levels", 10)

        contour_set = ax.contour(
            zi_contours,
            levels=levels,
            colors="white",
            linewidths=1
        )

        contour_file = os.path.join(
            working_dir, f"{map_choice}_contours.npz"
        )

        save_contours(
            contour_file,
            {f"level_{i}": c for i, c in enumerate(contour_set.allsegs)}
        )
    
    # Save figure
    os.makedirs(output_dir_path, exist_ok=True)
    output_path = os.path.join(
        output_dir_path, f"{map_choice}_map.pdf"
    )
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()

    print(f"Map saved in {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Plotting maps from analysed spectra using PyPistrello: Python Package for Integrating Spectral lines using Trapezoids, Error estimation and Line-features Optimization.')
    parser.add_argument('-t', '--input-file', type=str, required=True, help='FITS table with results from spectral lines analysis.')
    parser.add_argument('-c', '--config-file', type=str, required=True, help='Configuration YAML filename with parameters for plotting')
    parser.add_argument('-o', '--output-dir', type=str, required=True, help='Output directory to save results')
    parser.add_argument( "--map", type=str, required=True, choices=["flux", "vel", "snr"], help="Choose the type of map: flux, vel, snr" )
    args = parser.parse_args()

    fits_filename = args.input_file
    config_filename = args.config_file
    output_dir = args.output_dir
    map_choice = args.map
    
    print("\n")
    print(f"{BOLD}-----------------------------  PyPISTRELLO  ------------------------------")
    print("\U0001F987 Python Program for Integrating Spectral lines using TRapezoids,")
    print("Error estimation and Line-features Optimization \U0001F987 MAPS PLOTTING")
    print(f"--------------------------------------------------------------------------{RESET}")
    print("\n")

    working_dir = Path('.').resolve()
    print(f"Working directory: {working_dir}")

    # INPUT PROTECTIONS
    # protection against non-FITS files
    if not re.search(r'\.fits?$', fits_filename, re.IGNORECASE):
        raise ValueError(f"Input file '{fits_filename}' is not a FIT/FITS file. Please provide a valid FIT/FITS file.")
    
    fits_path = working_dir / fits_filename
    print(f"Input FITS file: {fits_path}")
    # make sure the input file exists
    if not fits_path.is_file():
        raise FileNotFoundError(f"Input file '{fits_filename}' does not exist. Please provide a valid file path.")
    
    # CONFIGURATION YAML FILE PROTECTIONS
    config_path = working_dir / config_filename
    print(f"Checking configuration file in {config_path}")
    if not config_path.is_file(): 
        raise FileNotFoundError(f"Configuration file '{config_filename}' does not exist. Please provide a valid file path." ) 
    if not re.search(r"\.ya?ml$", config_filename, re.IGNORECASE): 
        raise ValueError( f"Configuration file '{config_filename}' is not a YAML file. Please provide a valid YAML file." )
    
    # OUTPUT PROTECTIONS
    # make sure the output does not contain extensions because it is a directory
    if re.search(r'\.[a-zA-Z0-9]+$', output_dir): # this means there is a file extension: a "." followed by alphanumeric characters at the end of the string
        raise ValueError(f"Output directory '{output_dir}' should not contain file extensions.")
    
    # protection against overwriting existing output directory
    output_dir_path = working_dir / output_dir
    
    if output_dir_path.exists():
        print(f"WARNING: Output directory '{output_dir}' already exists.")
        response = input("Do you want to continue and overwrite existing files? (y/n): ")
        if response.lower() != "y": 
            print("Exiting program to prevent overwriting existing files.") 
            exit(0)
        else:
            print(f"{GREEN}INFO:{RESET} Using existing output directory: {output_dir}")
    else:
        output_dir_path.mkdir(parents=True) 
        print(f"{GREEN}INFO:{RESET} Created output directory: {output_dir}")

    print(f"{GREEN}INFO:{RESET} All inputs validated. Starting map generation...")
    map_plotting(working_dir, fits_path, config_path, output_dir_path, map_choice)


if __name__ == "__main__":
    main()