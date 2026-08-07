#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from astropy.io import fits
from pathlib import Path

import argparse
import matplotlib.pyplot as plt
import re

from .file_loading.load_yaml_file import load_yaml_file
from .file_loading.load_fits_table import load_fits_table
from .map_visualization.make_a_map import make_a_map

GREEN   = "\033[92m"
BLUE    = "\033[94m"
MAGENTA = "\033[95m"
BOLD = "\033[1m"
RESET   = "\033[0m"

def map_plotting(working_dir, fits_path, config_path, output_dir_path, map_choice): #bin_map_data=None
    """
    Generate and plot maps from spectral line analysis results.

    This function:
    1. Loads plotting and analysis parameters from a YAML configuration file.
    2. Reads a FITS table containing analysis results and its associated WCS.
    3. Calls the main plotting routine to generate the requested map.

    Parameters
    ----------
    working_dir : str
        Path to the working directory where intermediate files
        (e.g., contour files) will be saved.
    fits_path : str
        Path to the FITS file containing the analysis results table.
    config_path : str
        Path to the YAML configuration file with plotting parameters.
    output_dir_path : str
        Directory where the final map products will be stored.
    map_choice : str
        Identifier of the map to be generated (e.g., flux, velocity,
        dispersion). This value is used to select parameters and to
        name output files.

    Returns
    -------
    None
        This function does not return any value. It produces plots
        and saves output files to disk.
    """
    
    # Extract parameters from YAML configuration file:
    print(f"{BLUE}{BOLD} Reading parameters from YAML file{RESET}")
    config_parameters = load_yaml_file(config_path)
    print("YAML file read successfully")

    # Read the table with results from analysis
    table, wcs = load_fits_table(fits_path)

    # Plotting function
    make_a_map(table, wcs, config_parameters, working_dir, output_dir_path, map_choice) #bin_map_data
    print("INFO: Map created!")


def main():
    parser = argparse.ArgumentParser(description='Plotting maps from analysed spectra using PyPistrello: Python Package for Integrating Spectral lines using Trapezoids, Error estimation and Line-features Optimization.')
    parser.add_argument('-t', '--input-file', type=str, required=True, help='FITS table with results from spectral lines analysis.')
    parser.add_argument('-c', '--config-file', type=str, required=True, help='Configuration YAML filename with parameters for plotting')
    parser.add_argument('-o', '--output-dir', type=str, required=True, help='Output directory to save results')
    parser.add_argument( "--map", type=str, required=True, choices=["flux", "vel", "snr", "voronoi", "sigma", "EW"], help="Choose the type of map: flux, vel, snr, voronoi, sigma, EW" )
    args = parser.parse_args()

    fits_filename = args.input_file
    config_filename = args.config_file
    #bin_map = args.bin_map
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
    map_plotting(working_dir, fits_path, config_path, output_dir_path, map_choice) #bin_map_data


if __name__ == "__main__":
    main()