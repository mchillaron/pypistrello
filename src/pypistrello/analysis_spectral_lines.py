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

def analysis_spectral_lines(working_dir, fits_path, data_extension, output_dir_path, config_path, table_path,
                            simulation_dir_path, debug_level, run_voronoi):
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
    parser.add_argument('-vor', '--voronoi-tessellation', action='store_true', help='Whether to run Voronoi tessellation and binning of spectra with Powerbin.')
    parser.add_argument('-sim', '--simulations-dir', type=str, help='Directory name where simulated data is found.')
    parser.add_argument('--debug', type=int, default=0, help='Debug level (default: 0). Higher values (1 and 2) may print more detailed information for debugging purposes.')
    args = parser.parse_args()

    fits_filename = args.input_file
    data_extension = args.data_extension
    config_filename = args.config_file
    output_dir = args.output_dir
    table_file = args.table
    run_voronoi = args.voronoi_tessellation
    simulation_dir = args.simulations_dir
    debug_level = args.debug
    
    print("\n")
    #print(f"{BOLD}-----------------------------  PyPISTRELLO  ------------------------------")
    print(f"{BOLD} \U0001F987 Welcome to PyPISTRELLO \U0001F987")
    print(" Python Program for Integrating Spectral lines using TRapezoids,")
    print(f" Error estimation and Line-features Optimization {RESET}")
    #print(f"--------------------------------------------------------------------------{RESET}")
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
    print(f"Extracting data from extension {data_extension} of the FITS file.")
    
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
    print(f"The parameters measured will be saved in a FITS table named '{table_file}' in the output directory")
    table_path = output_dir_path / table_file

    # SIMULATIONS DIRECTORY
    # make sure the output does not contain extensions because it is a directory
    if simulation_dir is not None:
        if re.search(r'\.[a-zA-Z0-9]+$', simulation_dir):
            raise ValueError(f"Simulation directory '{simulation_dir}' should not contain file extensions.")
        simulation_dir_path = working_dir / simulation_dir
        # print a preview of the first few items found in the simulations directory
        if os.path.exists(simulation_dir_path):
            print(f"{GREEN}INFO:{RESET} Simulation directory '{simulation_dir}' found. Some of the files inside are:")
            files = sorted(os.listdir(simulation_dir_path))
            if len(files) <= 10:
                for item in files:
                    print(f"  - {item}")
            else:
                for item in files[:5]:
                    print(f"  - {item}")

                print("  ...")

                for item in files[-5:]:
                    print(f"  - {item}")
        else:
            raise ValueError(f"{GREEN}INFO:{RESET} Simulation directory '{simulation_dir}' not found. Please provide a valid directory path.")
    else:
        print(f"{GREEN}INFO:{RESET} No simulation directory provided. No simulation analysis will be performed for error estimation.")
        simulation_dir_path = None

    # VORONOI TESSELLATION
    if run_voronoi:
        print(f"{GREEN}INFO:{RESET} Voronoi tessellation and binning of spectra will be performed")
    else:
        print(f"{GREEN}INFO:{RESET} Voronoi tessellation and binning of spectra will not be performed")

    # DEBUG
    # DATA EXTENSION
    if debug_level < 0 and debug_level > 2:
        raise ValueError(f"Invalid debug number of '{debug_level}': must be 0, 1 or 2")
    print(f"Extracting data from extension {data_extension} of the FITS file.")


    analysis_spectral_lines(working_dir, fits_path, data_extension, output_dir_path, config_path, table_path, 
                            simulation_dir_path, debug_level, run_voronoi)


if __name__ == "__main__":
    main()