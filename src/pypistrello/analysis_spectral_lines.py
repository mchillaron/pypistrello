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

def analysis_spectral_lines(fits_path, output_dir_path):
    """Main function to analyze spectral lines from a FITS file and save results to an output directory.
    
    Parameters
    ----------
    fits_path : Path
        Path to the input FITS file containing spectral data.
    output_dir_path : Path
        Path to the output directory where results will be saved.
    """

    # Once we have the spectra in a FITS file and the directory in which we will be saving the results,
    # we can proceed to read the FITS file to extract the spectra.
    pass




def main():
    parser = argparse.ArgumentParser(description='Python Package for Integrating Spectral lines using Trapezoids, Error estimation and Line-features Optimization.')
    parser.add_argument('-F', '--input-file', type=str, required=True, help='Input file name containing spectral data, has to be FIT/FITS format')
    parser.add_argument('-o', '--output-dir', type=str, required=True, help='Output directory to save results')
    args = parser.parse_args()

    fits_name = args.input_file
    output_dir = args.output_dir

    # INPUT PROTECTIONS
    # protection against non-FITS files
    if not re.search(r'\.fits?$', fits_name, re.IGNORECASE):
        raise ValueError(f"Input file '{fits_name}' is not a FIT/FITS file. Please provide a valid FIT/FITS file.")
    
    fits_path = Path('.')/fits_name
    print(f"Input FITS file: {fits_path}")
    # make sure the input file exists
    if not fits_path.is_file():
        raise FileNotFoundError(f"Input file '{fits_name}' does not exist. Please provide a valid file path.")
    
    
    
    # OUTPUT PROTECTIONS
    # protection against overwriting existing output directory
    if os.path.exists(output_dir):
        raise FileExistsError(f"Output directory '{output_dir}' already exists. Please choose a different name to avoid overwriting existing data.")
    
    # make sure the output does not contain extensions becuase it is a directory
    if re.search(r'\.[a-zA-Z0-9]+$', output_dir): # this means there is a file extension: a "." followed by alphanumeric characters at the end of the string
        raise ValueError(f"Output directory '{output_dir}' should not contain file extensions.")

    output_dir_path = Path('.')/output_dir
    output_dir_path.mkdir(parents=True, exist_ok=True)
    print(f"Output files will be saved at: {output_dir_path}")

    analysis_spectral_lines(fits_path, output_dir_path)


if __name__ == "__main__":
    main()