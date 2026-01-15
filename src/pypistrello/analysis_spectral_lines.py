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

def analysis_spectral_lines(fits_path, 
                            output_dir_path, 
                            wavelength_path, 
                            config_path, 
                            redshift, 
                            line_restframe):
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
    """
 
    pass




def main():
    parser = argparse.ArgumentParser(description='Python Package for Integrating Spectral lines using Trapezoids, Error estimation and Line-features Optimization.')
    parser.add_argument('-F', '--input-file', type=str, required=True, help='Input file name containing coordinates and spectra, has to be FIT/FITS table format')
    parser.add_argument('-w', '--wavelength-range', type=str, required=True, help='Wavelength range to analyze, format: CSV')
    parser.add_argument('-c', '--config-file', type=str, required=True, help='Configuration YAML filename with parameters for analysis')
    parser.add_argument('-o', '--output-dir', type=str, required=True, help='Output directory to save results')
    parser.add_argument('-z', '--redshift', type=float, default=0.0, help='Redshift value to adjust spectral lines (default: 0.0)')
    parser.add_argument('-lrf', '--line-restframe', type=float, nargs='+', required=True, help='Rest-frame wavelength of spectral line to analyze')
    args = parser.parse_args()

    fits_filename = args.input_file
    wavelength_filename = args.wavelength_range
    config_filename = args.config_file
    output_dir = args.output_dir
    redshift = args.redshift
    line_restframe = args.line_restframe
    
    # INPUT PROTECTIONS
    # protection against non-FITS files
    if not re.search(r'\.fits?$', fits_filename, re.IGNORECASE):
        raise ValueError(f"Input file '{fits_filename}' is not a FIT/FITS file. Please provide a valid FIT/FITS file.")
    
    fits_path = Path('.')/fits_filename
    print(f"Input FITS file: {fits_path}")
    # make sure the input file exists
    if not fits_path.is_file():
        raise FileNotFoundError(f"Input file '{fits_filename}' does not exist. Please provide a valid file path.")
    
    # OUTPUT PROTECTIONS
    # protection against overwriting existing output directory
    if os.path.exists(output_dir):
        raise FileExistsError(f"Output directory '{output_dir}' already exists. Please choose a different name to avoid overwriting existing data.")
    
    # make sure the output does not contain extensions because it is a directory
    if re.search(r'\.[a-zA-Z0-9]+$', output_dir): # this means there is a file extension: a "." followed by alphanumeric characters at the end of the string
        raise ValueError(f"Output directory '{output_dir}' should not contain file extensions.")

    output_dir_path = Path('.')/output_dir
    output_dir_path.mkdir(parents=True, exist_ok=True)
    print(f"Output files will be saved at: {output_dir_path}")

    wavelength_path = Path('.')/wavelength_filename
    print(f"Wavelength range file: {wavelength_path}")

    config_path = Path('.')/config_filename
    print(f"Configuration file: {config_path}")

    analysis_spectral_lines(fits_path, output_dir_path, wavelength_path, config_path, redshift, line_restframe)


if __name__ == "__main__":
    main()