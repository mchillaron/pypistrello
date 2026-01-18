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

from .load_fits_table import load_fits_table
from .load_wavelength_range import load_wavelength_range
from .diagnostic_plot.plot_diagnostic_spectra import plot_diagnostic_spectra

def analysis_spectral_lines(fits_path, 
                            output_dir_path, 
                            wavelength_path, 
                            config_path, 
                            redshift, 
                            line_restframe, 
                            diagnostic_spectra):
    """Main function to analyze spectral lines from a FITS file and save results to an output directory.
    
    Parameters
    ----------
    fits_path : Path
        Path to the input FITS file containing a table with coordinates and spectra.
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
    
    # load the FITS table from fits_path
    spectra_table = load_fits_table(fits_path)
    print(f"INFO: FITS table loaded from {fits_path}. Check below for a preview:")
    print(spectra_table[:10])  # show in the terminal the first 20 rows of the table
    # show also the header of the table
    # print(spectra_table.meta)

    # load the wavelength range from wavelength_path
    wavelength_range = load_wavelength_range(wavelength_path)
    print(f"INFO: Wavelength range loaded from {wavelength_path}.")

    # Diagnostic plot to prepare analysis
    if diagnostic_spectra is not None:
        plot_diagnostic_spectra(spectra_table, wavelength_range, diagnostic_spectra, output_dir_path, redshift, line_restframe)
    





def main():
    parser = argparse.ArgumentParser(description='Python Package for Integrating Spectral lines using Trapezoids, Error estimation and Line-features Optimization.')
    parser.add_argument('-F', '--input-file', type=str, required=True, help='Table filename containing coordinates and spectra, has to be FIT/FITS table format')
    parser.add_argument('-w', '--wavelength-range', type=str, required=True, help='Wavelength range to analyze, format: CSV')
    parser.add_argument('-c', '--config-file', type=str, required=True, help='Configuration YAML filename with parameters for analysis')
    parser.add_argument('-o', '--output-dir', type=str, required=True, help='Output directory to save results')
    parser.add_argument('-d', '--diagnostic-spectra', type=int, nargs=4, metavar=("x1", "x2", "y1", "y2"), help="Coordinates of spectra to integrate for diagnostic plot: x1 x2 y1 y2. FITS indices.")
    parser.add_argument('-z', '--redshift', type=float, default=0.0, help='Redshift value to adjust spectral lines (default: 0.0)')
    parser.add_argument('-lrf', '--line-restframe', type=float, nargs='+', required=True, help='Rest-frame wavelength of spectral line to analyze')
    args = parser.parse_args()

    fits_filename = args.input_file
    wavelength_filename = args.wavelength_range
    config_filename = args.config_file
    output_dir = args.output_dir
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
    
    # WAVELENGTH FILE
    wavelength_path = working_dir / wavelength_filename
    if not os.path.isfile(wavelength_filename):
        raise FileNotFoundError(f"Wavelength range file '{wavelength_filename}' does not exist. Please provide a valid file path.")
    if not re.search(r'\.csv$', wavelength_filename, re.IGNORECASE):
        raise ValueError(f"Wavelength range file '{wavelength_filename}' is not a CSV file. Please provide a valid CSV file.")
    
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
    print(type(line_restframe))
    print(line_restframe)

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
        if x2 <= x1 or y2 <= y1:
            raise ValueError(f"Diagnostic spectra coordinates are invalid. Ensure that x2 > x1 and y2 > y1.")
        
        print(f"For the diagnostic plot, the spectra will be integrated over the provided coordinates: {x1, x2, y1, y2}")
    

    analysis_spectral_lines(fits_path, output_dir_path, wavelength_path, config_path, redshift, line_restframe, diagnostic_spectra)


if __name__ == "__main__":
    main()