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

import numpy as np
import os
import re
import sys
import time

from .file_loading.load_fits_cube import read_fits_cube
from .file_loading.load_wavelength_range import load_wavelength_range
from .file_loading.load_yaml_file import load_yaml_file
from .file_loading.get_wavelength_axis import get_wavelength_axis
from .file_loading.save_table_fits import save_table_with_wcs_extension
from .file_loading.yn_question import question_yes_no

from .diagnostic_plot.plot_diagnostic_spectra import plot_diagnostic_spectra
from .area_fitting.main_trapz_fitting import main_trapz_fitting
from .area_fitting.area_trapz_spectra_bin import area_trapz_spectra_bin
from .analysis_tools.measure_spectra_properties import measure_spectra_properties

from .simulated_data.process_simulations import process_simulations
from .simulated_data.compute_sim_snr import compute_sim_snr
from .voronoi_binning.run_powerbin import run_powerbin
from .voronoi_binning.sum_spectra_voronoi import sum_spectra_voronoi
from .voronoi_binning.build_voronoi_table import build_voronoi_table
from .voronoi_binning.extract_spectra_from_table import extract_spectra_from_table
from .voronoi_binning.extract_spectra_from_table import check_alignment
from .voronoi_binning.propagate_bin_to_spaxel import propagate_bin_to_spaxel_table
from .voronoi_binning.save_voronoi_bin_spectra import save_voronoi_bin_spectra_pdf


GREEN   = "\033[92m"
BLUE    = "\033[94m"
MAGENTA = "\033[95m"
BOLD = "\033[1m"
RESET   = "\033[0m"

def analysis_spectral_lines(working_dir, fits_path, data_extension, output_dir_path, config_path, table_path,
                            simulation_dir_path, debug_level, run_voronoi, start):
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
    
    # load the FITS datacube and information from headers
    print(f"{BLUE}{BOLD} Reading header and data from FITS cube{RESET}")
    primary_header_real, data_header_real, cube_data_real, wcs_info_real = read_fits_cube(fits_path, data_extension)
    print("Cube headers and data read successfully")

    # Load de YAML file and read parameters
    print(f"{BLUE}{BOLD} Reading parameters from YAML file{RESET}")
    config_parameters = load_yaml_file(config_path)
    print("YAML file read successfully")

    line_name = config_parameters["line_name"]
    print(f"Analysing {line_name} line")

    line_restframe = config_parameters["line_restframe"]
    if not all(isinstance(lrf, float) for lrf in line_restframe):
        raise ValueError(f"One or more line rest-frame wavelengths are not floats. Please provide valid float values using [].")
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
        wavelength_range = get_wavelength_axis(data_header_real, primary_header_real)
        print(f"{GREEN}INFO:{RESET} Wavelength range calculated.")

    # Ask whether to run the diagnostic plot or not
    run_diag = question_yes_no("Do you want to run the diagnostic diagram?")

    if run_diag:
        print("INFO: Running diagnostic diagram...")
        plot_diagnostic_spectra(cube_data_real, wavelength_range, output_dir_path, config_parameters, redshift, line_restframe)
        print("INFO: Diagnostic diagram completed.")
        sys.exit(0)

    else:
        print("INFO: Skipping diagnostic diagram.")
        print(f"{BLUE}{BOLD} Integrating {line_name} line area with trapezoids {RESET}")

        real_cube_measured = False              # creat a flag that will be useful when analysing simulated cubes

        table_results_fitting = main_trapz_fitting(output_dir_path, cube_data_real, wcs_info_real,
                                                  wavelength_range, config_parameters, table_path,
                                                  redshift, line_restframe)
        
        table_results_fitting[:5].pprint()
        table_results_fitting[-5:].pprint()
        
        if simulation_dir_path is not None:
            print(f"{BLUE}{BOLD} Working with simulated cubes{RESET}")
            sim_file = output_dir_path / "simulated_measurements.npy"
            sim_file_npz = output_dir_path / "simulated_measurements.npz"

            if sim_file.exists():
                load_simulated_data = question_yes_no(
                    "Load simulated measurements from existing .npy file?"
                )
                if load_simulated_data:
                    print(f"{GREEN}INFO:{RESET} Loading simulated measurements from {sim_file}")
                    simulation_results = np.load(sim_file)
                else:
                    print(f"{GREEN}INFO:{RESET} Rerunning measurements on simulated cubes in {simulation_dir_path}")
                    simulation_results = process_simulations(simulation_dir_path, output_dir_path, wavelength_range,
                                                            data_extension, config_parameters, redshift, line_restframe)
            else:
                print(f"{GREEN}INFO:{RESET} Looking for simulated data in {simulation_dir_path}")
                simulation_results = process_simulations(simulation_dir_path, output_dir_path, wavelength_range,
                                                        data_extension, config_parameters, redshift, line_restframe)

            print(f"{BLUE}{BOLD} Calculating SNR using simulated measurements{RESET}")
            snr_table = compute_sim_snr(table_results_fitting, "area_trapz", simulation_results, config_parameters, debug_level=debug_level)
            table_results_fitting["snr"] = snr_table

        else:
            print(f"{GREEN}INFO:{RESET} No simulation directory provided. Skipping SNR calculation from simulations.")
            snr_table = None
            print(f"{MAGENTA} Attention!{RESET} The SNR values have been calculated as the ratio of the line flux ")
            print(f" to the noise in the continuum, both estimated from the real data without using simulations.")
            print(f" These SNR values are less accurate than those calculated using simulations and UNDERESTIMATE the true SNR,")
            print(f" they only provide a rough estimate of the SNR for each spectrum.")


        if run_voronoi:
            print(f"{BLUE}{BOLD} Applying Powerbin for Voronoi Tessellation {RESET}")

            pow, table_results_fitting, pow_valid_mask = run_powerbin(table_results_fitting, config_parameters, debug_level, snr_table=snr_table)

            cube_binned, bin_map, cube_voronoi = sum_spectra_voronoi(cube_data_real, table_results_fitting, output_dir_path, debug_level)
            hdu = fits.PrimaryHDU(cube_binned)
            hdu.writeto(output_dir_path / "cube_binned.fits", overwrite=True)

            # adapt the Table to the summed spectra after voronoi binning
            analysis_table = build_voronoi_table(table_results_fitting, pow)
            spectra = cube_binned   # (n_lambda, n_bins)

            if debug_level > 1:
                save_voronoi_bin_spectra_pdf(
                    spectra=cube_binned,
                    wavelength=wavelength_range,
                    analysis_table=analysis_table,
                    config_parameters=config_parameters,
                    output_dir=output_dir_path,
                    simulated_spectra=None,  
                    grid_size=(5, 5),
                    sort_by_snr=True
                )
            
            analysis_table = area_trapz_spectra_bin(spectra, wavelength_range, config_parameters, 
                                                analysis_table, redshift, line_restframe, debug_level)
            
        else:
            print(f"{GREEN}INFO:{RESET} No Voronoi binning")

            spectra = extract_spectra_from_table(cube_data_real, table_results_fitting)
            analysis_table = table_results_fitting
            cube_binned = None
            bin_map = None
            cube_voronoi = None

            if check_alignment(cube_data_real, analysis_table, spectra):
                print(f"{GREEN}INFO:{RESET} Spectra alignment check PASSED")
            else:
                print(f"{MAGENTA}WARNING:{RESET} Spectra alignment check FAILED. There may be a mismatch between the spectra extracted from the cube and the coordinates in the table.")
                ValueError("Please check the coordinates in the table and the structure of the cube data.")


        print("Data preparation ready for analysis of spectral lines") 
        analysis_table = measure_spectra_properties(spectra, wavelength_range, config_parameters,
                                                    analysis_table, redshift, line_restframe, output_dir_path, real_cube_measured, debug_level)
        real_cube_measured = True           # change the flag once the real cube has been analysed
        
        # Saving the analysis table with information for every spaxel
        if run_voronoi:
            columns_to_copy = ["n_pix", "bin_center_x", "bin_center_y", 
                               "bin_area_trapz","bin_cont_noise","bin_snr_trapz","bin_cont_coeffs",  #"bin_snr_simulated",
                               "velocity", "offsets",
                               "amp_gauss", "mu_gauss", "sigma_gauss", "fwhm", "cont_gauss", "area_gauss", "chi2_gauss",
                               "amp_ha", "amp_nii6548", "amp_nii6583", "area_ha", "area_nii6548", "area_nii6583", "area_total,"
                               "amp1", "amp2", "mu1", "mu2", "area1", "area2", "area_total"]
            
            table_collapsed = propagate_bin_to_spaxel_table(table_results_fitting, analysis_table, columns_to_copy)
        else:
            table_collapsed = analysis_table
        
        save_table_with_wcs_extension(
            table_collapsed,
            table_path,
            wcs_info=wcs_info_real
        )
        table_collapsed[:5].pprint()
        table_collapsed[-5:].pprint()
        
        print("\n" + "=" * 70)
        print("             TABLE SUMMARY")
        print("=" * 70)

        print(f"Number of rows    : {len(table_collapsed):>8}")
        print(f"Number of columns : {len(table_collapsed.colnames):>8}")

        print("\nColumn names")
        print("-" * 70)

        for i, col in enumerate(table_collapsed.colnames, start=1):
            print(f"{i:2d}. {col}")

        print("=" * 70)
        print("INFO: Final unified table saved")

        # -------------------------
        # If simulations are provided, the same parameters are measured in every cube and
        # statistics are saved in final file.

        if simulation_dir_path is not None:
            print(f"{BLUE}{BOLD} Performing same analysis on simulated cubes {RESET}")
            if sim_file_npz.exists():
                print("A .npz file with metadata has been found in the directory")
                data_sim_npz = np.load(sim_file_npz, allow_pickle=True)

                if "columns" in data_sim_npz:
                    columns = data_sim_npz["columns"].tolist()
                    print(f"{BLUE}Columns in simulation file:{RESET}")
                    for col in columns:
                        print(f"  - {col}")
                else:
                    print(f"{MAGENTA}WARNING:{RESET} No column metadata found in .npz")
                
                load_simulated_data = question_yes_no(
                    "Load simulated measurements from existing .npy file?"
                )
                if load_simulated_data:
                    print(f"{GREEN}INFO:{RESET} Loading simulated measurements from {sim_file}")
                    simulation_results_props = np.load(sim_file)
                else:
                    print(f"{GREEN}INFO:{RESET} Rerunning measurements on simulated cubes in {simulation_dir_path}")
                    simulation_results_props = process_simulations(simulation_dir_path, output_dir_path, wavelength_range,
                                                        data_extension, config_parameters, redshift, line_restframe,
                                                        real_cube_measured, snr_table=snr_table, pow=pow, pow_valid_mask=pow_valid_mask)
            else:
                print("INFO: Creating a new .npy file with all measurements from simulated cubes...")
                simulation_results_props = process_simulations(simulation_dir_path, output_dir_path, wavelength_range,
                                                            data_extension, config_parameters, redshift, line_restframe,
                                                            real_cube_measured, snr_table=snr_table, pow=pow, pow_valid_mask=pow_valid_mask)

            if run_voronoi:
                # Adding to every table a new column called "bin_snr_sim" 
                bin_snr_table = compute_sim_snr(table_collapsed, "bin_area_trapz", simulation_results_props, config_parameters, debug_level=debug_level)
                if len(table_collapsed) == len(bin_snr_table):
                    col_index = table_collapsed.colnames.index("bin_snr_trapz") + 1 
                    table_collapsed.add_column(bin_snr_table, name="bin_snr_sim", index=col_index)
                    table_collapsed.write(table_path, overwrite=True)
        
                
    print(f"Goodbye!")
    end = time.perf_counter()
    elapsed = end - start
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    seconds = elapsed % 60
    print(f"⏱️ Execution time: {elapsed:.3f} s ({hours:02d}:{minutes:02d}:{seconds:06.3f})")
    

            


def main():
    start = time.perf_counter()
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
    if debug_level < 0 or debug_level > 2:
        raise ValueError(f"Invalid debug number of '{debug_level}': must be 0, 1 or 2")
    print(f"Extracting data from extension {data_extension} of the FITS file.")


    analysis_spectral_lines(working_dir, fits_path, data_extension, output_dir_path, config_path, table_path, 
                            simulation_dir_path, debug_level, run_voronoi, start)


if __name__ == "__main__":
    main()