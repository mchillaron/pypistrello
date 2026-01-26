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
import os
import re

from .file_loading.load_fits_table import load_fits_table
from .file_loading.load_fits_cube import read_fits_cube
from .file_loading.load_wavelength_range import load_wavelength_range
from .file_loading.load_yaml_file import load_yaml_file
from .file_loading.get_wavelength_axis import get_wavelength_axis
from .file_loading.load_yaml_file import validate_region_config
from .diagnostic_plot.plot_diagnostic_spectra import plot_diagnostic_spectra
from .line_fitting.main_line_fitting import main_line_fitting

def analysis_spectral_lines(fits_path,
                            data_extension,
                            output_dir_path, 
                            wavelength_path, 
                            config_path, 
                            table_parameters_path,
                            redshift, 
                            line_restframe, 
                            diagnostic_spectra):
    """Main function to analyze spectral lines from a FITS file and save results to an output directory.
    
    Parameters
    ----------
    fits_path : Path
        Path to the input FITS file containing a table with coordinates and spectra.
    data_extension: int
        Extension number of the FITS cube where data is found.
    output_dir_path : Path
        Path to the output directory where results will be saved.
    wavelength_path : Path
        Path to the file containing the wavelength range to analyze.
    config_path : Path
        Path to the configuration YAML file with parameters for analysis.
    redshift : float
        Redshift value to adjust spectral lines.
    line_restframe : list of float
        List of rest-frame wavelengths of spectral lines to analyze.
    diagnostic_spectra : tuple of int or None
        Coordinates of spectra to integrate for diagnostic plot (x1, x2, y1, y2) or None if no diagnostic plot is needed.
    """
    
    if data_extension == 0:
        print("INFO: Using extension 0 as data_header")
    else:
        print(f"INFO: Using extension 0 as primary_header and extension {data_extension} as data_header")

    # load the FITS datacube and information from headers
    primary_header, data_header, cube_data = read_fits_cube(fits_path, data_extension)

    # load the wavelength range from wavelength_path
    if wavelength_path is not None:
        wavelength_range = load_wavelength_range(wavelength_path)
        print(f"INFO: Wavelength range loaded from {wavelength_path}.")
    else:
        wavelength_range = get_wavelength_axis(data_header, primary_header)
        print("INFO: Wavelength range calculated")
        print(wavelength_range)

    # Diagnostic plot to prepare analysis
    if diagnostic_spectra is not None:
        plot_diagnostic_spectra(cube_data, wavelength_range, diagnostic_spectra, output_dir_path, redshift, line_restframe)
    
    # After the diagnostic plot, we already have all the parameters to fill in the YAML file.
    config_parameters = load_yaml_file(config_path)

    tab_line_fit = main_line_fitting(output_dir_path, cube_data, wavelength_range, config_parameters,
                     table_parameters_path, redshift, line_restframe)
    print(tab_line_fit)
    





def main():
    parser = argparse.ArgumentParser(description='Python Package for Integrating Spectral lines using Trapezoids, Error estimation and Line-features Optimization.')
    parser.add_argument('-F', '--input-file', type=str, required=True, help='FITS datacube filename to work with.')
    parser.add_argument('-ext', '--data-extension', type=int, required=True, help='FITS extension number where data is found (>= 0).')
    parser.add_argument('-w', '--wavelength-range', type=str, help='Wavelength range to analyze, format: CSV')
    parser.add_argument('-c', '--config-file', type=str, required=True, help='Configuration YAML filename with parameters for analysis')
    parser.add_argument('-o', '--output-dir', type=str, required=True, help='Output directory to save results')
    parser.add_argument('-t', '--table-parameters', type=str, help='Name of FITS table to save the line-fit parameters. Include .fits extension.')
    parser.add_argument('-d', '--diagnostic-spectra', type=int, nargs=4, metavar=("x1", "x2", "y1", "y2"), help="Coordinates of spectra to integrate for diagnostic plot: x1 x2 y1 y2. FITS indices.")
    parser.add_argument('-z', '--redshift', type=float, required=True, default=0.0, help='Redshift value to adjust spectral lines (default: 0.0)')
    parser.add_argument('-lrf', '--line-restframe', type=float, nargs='+', required=True, help='Rest-frame wavelength of spectral line to analyze')
    args = parser.parse_args()

    fits_filename = args.input_file
    data_extension = args.data_extension
    wavelength_filename = args.wavelength_range
    config_filename = args.config_file
    output_dir = args.output_dir
    table_parameters = args.table_parameters
    diagnostic_spectra = args.diagnostic_spectra
    redshift = args.redshift
    line_restframe = args.line_restframe
    
    print("\n")
    print("-----------------------------  PyPISTRELLO  ------------------------------")
    print("\U0001F987 Python Program for Integrating Spectral lines using TRapezoids,")
    print("Error estimation and Line-features Optimization \U0001F987")
    print("--------------------------------------------------------------------------")
    print("\n")


    working_dir = Path('.').resolve()
    print(f"Working directory: {working_dir}")
    print(f"INFO: All files will be read/written relative to the working directory.")


    # INPUT PROTECTIONS
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
    
    # OUTPUT PROTECTIONS
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
            print(f"INFO: Using existing output directory: {output_dir}")
    else:
        os.makedirs(output_dir)
        print(f"INFO: Created output directory: {output_dir}")
    
    # CONFIGURATION YAML FILE PROTECTIONS
    config_path = working_dir / config_filename
    if not os.path.isfile(config_filename):
        raise FileNotFoundError(f"Configuration file '{config_filename}' does not exist. Please provide a valid file path.")
    if not re.search(r'\.ya?ml$', config_filename, re.IGNORECASE):
        raise ValueError(f"Configuration file '{config_filename}' is not a YAML file. Please provide a valid YAML file.")
    
    # TABLE PARAMETERS
    # protection against non-FITS files
    if not re.search(r'\.fits?$', table_parameters, re.IGNORECASE):
        raise ValueError(f"Input file '{table_parameters}' is not a FIT/FITS file. Please provide a valid FIT/FITS file.")
    
    table_parameters_path = output_dir_path / table_parameters
    
    # WAVELENGTH FILE
    if wavelength_filename is not None:
        wavelength_path = working_dir / wavelength_filename
        if not os.path.isfile(wavelength_filename):
            raise FileNotFoundError(f"Wavelength range file '{wavelength_filename}' does not exist. Please provide a valid file path.")
        if not re.search(r'\.csv$', wavelength_filename, re.IGNORECASE):
            raise ValueError(f"Wavelength range file '{wavelength_filename}' is not a CSV file. Please provide a valid CSV file.")
    else:
        wavelength_path = None
        
    # REDHSIFT PROTECTIONS
    if not isinstance(redshift, float):
        raise ValueError(f"Redshift value '{redshift}' is not a float. Please provide a valid float value.")
    if redshift < 0.0:
        raise ValueError(f"Redshift value '{redshift}' is negative. Please provide a non-negative float value.")
    if redshift == 0.0:
        print(f"INFO: Redshift value is 0.0, no adjustment will be made to spectral lines.")
    print(f"Redshift value: {redshift}")
    
    # LINE REST-FRAME PROTECTIONS
    if not all(isinstance(lrf, float) for lrf in line_restframe):
        raise ValueError(f"One or more line rest-frame wavelengths are not floats. Please provide valid float values.")
    print(f"Line rest-frame wavelengths: {line_restframe}")

    # DIAGNOSTIC SPECTRA PROTECTIONS
    if diagnostic_spectra is None:
        print(f"INFO: No diagnostic spectra coordinates provided, no diagnostic plot will be generated.")
    else:
        x1, x2, y1, y2 = diagnostic_spectra
        if len(diagnostic_spectra) != 4:
            raise ValueError(f"Diagnostic spectra coordinates must contain exactly 4 integers for coordinates of spectra to be plotted: x1, x2, y1, y2.")
        if not all(isinstance(coord, int) for coord in diagnostic_spectra):
            raise ValueError(f"One or more diagnostic spectra coordinates are not integers. Please provide valid integer values.")
        print(f"Diagnostic spectra coordinates: {diagnostic_spectra}")

        if x1 < 0 or x2 < 0 or y1 < 0 or y2 < 0:
            raise ValueError(f"Diagnostic spectra coordinates must be non-negative integers.")
        if x2 < x1 or y2 < y1:
            raise ValueError(f"Diagnostic spectra coordinates are invalid. Ensure that x2 >= x1 and y2 >= y1.")
        
        print(f"For the diagnostic plot, the spectra will be integrated over the provided coordinates: {x1, x2, y1, y2}")
    

    analysis_spectral_lines(fits_path, data_extension, output_dir_path, 
                            wavelength_path, config_path, table_parameters_path,
                            redshift, line_restframe, diagnostic_spectra)


if __name__ == "__main__":
    main()