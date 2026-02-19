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
import astropy.units as u
import os
import re
import sys

from .file_loading.load_fits_cube import read_fits_cube
from .file_loading.load_wavelength_range import load_wavelength_range
from .file_loading.load_yaml_file import load_yaml_file
from .file_loading.get_wavelength_axis import get_wavelength_axis
from .file_loading.save_table_fits import save_table_with_wcs_extension
from .file_loading.yn_question import question_yes_no

from .diagnostic_plot.plot_diagnostic_spectra import plot_diagnostic_spectra

from .line_fitting.crosscorrelation_spectra import crosscorrelate_spectra
from .line_fitting.run_powerbin import run_powerbin
from .line_fitting.sum_spectra_voronoi import sum_spectra_voronoi
from .line_fitting.crosscorrelation_spectra import convert_offset_velocity
from .line_fitting.main_line_fitting import main_line_fitting


GREEN   = "\033[92m"
BLUE    = "\033[94m"
MAGENTA = "\033[95m"
BOLD = "\033[1m"
RESET   = "\033[0m"

def analysis_spectral_lines(working_dir, fits_path, data_extension, output_dir_path, config_path, table_path):
    """Main function to analyze spectral lines from a FITS file and save results to an output directory.
    
    Parameters
    ----------
    fits_path : Path
        Path to the input FITS file containing a table with coordinates and spectra.
    data_extension: int
        Extension number of the FITS cube where data is found.
    output_dir_path : Path
        Path to the output directory where results will be saved.
    config_path : Path
        Path to the configuration YAML file with parameters for analysis.
    """
    
    if data_extension == 0:
        print(f"{GREEN}INFO:{RESET} Using extension 0 as data_header")
    else:
        print(f"{GREEN}INFO:{RESET} Using extension 0 as primary_header and extension {data_extension} as data_header")

    # load the FITS datacube and information from headers
    print(f"{BLUE}{BOLD} Reading header and data from FITS cube{RESET}")
    primary_header, data_header, cube_data, wcs_info = read_fits_cube(fits_path, data_extension)
    print("Cube headers and data read successfully")

    # Load de YAML file and read parameters
    print(f"{BLUE}{BOLD} Reading parameters from YAML file{RESET}")
    config_parameters = load_yaml_file(config_path)
    print("YAML file read successfully")

    line_name = config_parameters["line_name"]
    print(f"Analysing {line_name} line")

    line_restframe = config_parameters["line_restframe"]
    if not all(isinstance(lrf, float) for lrf in line_restframe):
        raise ValueError(f"One or more line rest-frame wavelengths are not floats. Please provide valid float values.")
    print(f"Line rest-frame wavelength: {line_restframe}")

    redshift = config_parameters["redshift"]
    if not isinstance(redshift, float):
        raise ValueError(f"Redshift value '{redshift}' is not a float. Please provide a valid float value.")
    print(f"Redshift value provided: {redshift}")

    # load the wavelength range from wavelength_path
    wavelength_param = config_parameters.get("wavelength_file")
    if wavelength_param is not None:
        wavelength_path = Path(working_dir) / wavelength_param
        if not os.path.isfile(wavelength_path):
            raise FileNotFoundError(f"Wavelength range file '{wavelength_path}' does not exist. Please provide a valid file path.")
        if not re.search(r'\.csv$', wavelength_param, re.IGNORECASE):
            raise ValueError(f"Wavelength range file '{wavelength_param}' is not a CSV file. Please provide a valid CSV file.")

        wavelength_range = load_wavelength_range(wavelength_path)
        print(f"{GREEN}INFO:{RESET} Wavelength range loaded from {wavelength_path}.")
    else:
        wavelength_range = get_wavelength_axis(data_header, primary_header)
        print(f"{GREEN}INFO:{RESET} Wavelength range calculated.")

    # Ask whether to run the diagnostic plot or not
    run_diag = question_yes_no("Do you want to run the diagnostic diagram?")

    if run_diag:
        print("INFO: Running diagnostic diagram...")
        plot_diagnostic_spectra(cube_data, wavelength_range, output_dir_path, config_parameters, redshift, line_restframe)
        print("INFO: Diagnostic diagram completed.")
        sys.exit(0)

    else:
        print("INFO: Skipping diagnostic diagram.")

        print(f"{BLUE}{BOLD} Integrating {line_name} line area with trapezoids {RESET}")
        table_results_fitting = main_line_fitting(output_dir_path, cube_data, wcs_info,
                                                  wavelength_range, config_parameters, table_path,
                                                  redshift, line_restframe)

        print(f"{BLUE}{BOLD} Applying Powerbin for Voronoi Tessellation {RESET}")
        pow, table_results_fitting = run_powerbin(table_results_fitting, config_parameters)

        print("Summing spectra in each Voronoi bin and creating bin_map and cube_voronoi")
        cube_2d_binned, bin_map, cube_voronoi = sum_spectra_voronoi(cube_data, table_results_fitting)

        # At this point, the calculations are carried out on Voronoi-binned spectra
        print(f"{BLUE}{BOLD} Crosscorrelation to reference spectrum {RESET}")
        offsets_pixel_array, fpeaf_croscorr_array = crosscorrelate_spectra(cube_2d_binned, wavelength_range, 
                                                                           config_parameters, table_results_fitting, 
                                                                           bin_map=bin_map)

        print(f"{BLUE}{BOLD} Calculating velocities for line {line_name} {RESET}")
        velocity = convert_offset_velocity(offsets_pixel_array, wavelength_range, redshift, line_restframe)

        table_results_fitting["offsets"] = offsets_pixel_array
        table_results_fitting["velocity"] = velocity * u.km / u.s
        table_results_fitting[:5].pprint()
        table_results_fitting[-5:].pprint()

        print(f"Saving velocity values to {table_path}")
        save_table_with_wcs_extension(table_results_fitting, table_path, wcs_info=wcs_info)
    
    


def main():
    parser = argparse.ArgumentParser(description='Python Package for Integrating Spectral lines using Trapezoids, Error estimation and Line-features Optimization.')
    parser.add_argument('-F', '--input-file', type=str, required=True, help='FITS datacube filename to work with.')
    parser.add_argument('-ext', '--data-extension', type=int, required=True, help='FITS extension number where data is found (>= 0).')
    parser.add_argument('-c', '--config-file', type=str, required=True, help='Configuration YAML filename with parameters for analysis')
    parser.add_argument('-o', '--output-dir', type=str, required=True, help='Output directory to save results')
    parser.add_argument('-t', '--table', type=str, help='Name of FITS table where to save the line-fit parameters. Include .fits extension.')
    parser.add_argument('--simulate', action='store_true', help='Run the program with simulated data instead of real FITS input.')
    args = parser.parse_args()

    fits_filename = args.input_file
    data_extension = args.data_extension
    config_filename = args.config_file
    output_dir = args.output_dir
    table_file = args.table
    simulate = args.simulate
    
    print("\n")
    print(f"{BOLD}-----------------------------  PyPISTRELLO  ------------------------------")
    print("\U0001F987 Python Program for Integrating Spectral lines using TRapezoids,")
    print("Error estimation and Line-features Optimization \U0001F987")
    print(f"--------------------------------------------------------------------------{RESET}")
    print("\n")


    working_dir = Path('.').resolve()
    print(f"Working directory: {working_dir}")
    print(f"{GREEN}INFO:{RESET} All files will be read/written relative to the working directory.")


    # INPUT
    # protection against non-FITS files
    if not re.search(r'\.fits?$', fits_filename, re.IGNORECASE):
        raise ValueError(f"Input file '{fits_filename}' is not a FIT/FITS file. Please provide a valid FIT/FITS file.")
    
    fits_path = working_dir / fits_filename
    print(f"Input FITS file: {fits_path}")
    # make sure the input file exists
    if not fits_path.is_file():
        raise FileNotFoundError(f"Input file '{fits_filename}' does not exist. Please provide a valid file path.")
    
    # DATA EXTENSION
    if data_extension < 0:
        raise ValueError(f"Invalid extension number of '{data_extension}': must be >= 0")
    
    # OUTPUT
    # make sure the output does not contain extensions because it is a directory
    if re.search(r'\.[a-zA-Z0-9]+$', output_dir): # this means there is a file extension: a "." followed by alphanumeric characters at the end of the string
        raise ValueError(f"Output directory '{output_dir}' should not contain file extensions.")
    
    # protection against overwriting existing output directory
    output_dir_path = working_dir / output_dir
    if os.path.exists(output_dir):
        print(f"WARNING: Output directory '{output_dir}' already exists.")                              # ask user for confirmation to overwrite
        response = input("Do you want to continue and overwrite existing files? (y/n): ")
        if response.lower() != 'y':
            print("Exiting program to prevent overwriting existing files.")
            exit(0)
        else:
            print(f"{GREEN}INFO:{RESET} Using existing output directory: {output_dir}")
    else:
        os.makedirs(output_dir)
        print(f"{GREEN}INFO:{RESET} Created output directory: {output_dir}")
    
    # CONFIGURATION YAML FILE
    config_path = working_dir / config_filename
    print(f"Checking configuration file in {config_path}")
    if not os.path.isfile(config_filename):
        raise FileNotFoundError(f"Configuration file '{config_filename}' does not exist. Please provide a valid file path.")
    if not re.search(r'\.ya?ml$', config_filename, re.IGNORECASE):
        raise ValueError(f"Configuration file '{config_filename}' is not a YAML file. Please provide a valid YAML file.")
    
    # TABLE PARAMETERS
    # protection against non-FITS files
    if not re.search(r'\.fits?$', table_file, re.IGNORECASE):
        raise ValueError(f"Input file '{table_file}' is not a FIT/FITS file. Please provide a valid FIT/FITS file.")

    table_path = output_dir_path / table_file

    analysis_spectral_lines(working_dir, fits_path, data_extension, output_dir_path, config_path, table_path)


if __name__ == "__main__":
    main()