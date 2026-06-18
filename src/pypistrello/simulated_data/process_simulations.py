#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from pathlib import Path
import numpy as np

from .process_simulated_cube import process_simulated_cube

GREEN   = "\033[92m"
RESET   = "\033[0m"
BOLD = "\033[1m"

def process_simulations(
        simulation_dir,
        output_dir_path,
        wavelength_range,
        data_extension,
        config_parameters,
        redshift,
        line_restframe,
        real_cube_measured=False, 
        snr_table=None,
        pow=None):

    cube_files = sorted(
        f for f in simulation_dir.iterdir()
        if f.suffix.lower() in [".fits", ".fit"]
    )
    print(f"{GREEN}INFO:{RESET} A total number of {len(cube_files)} simulation cubes will be processed")
    if len(cube_files) == 0:
        ValueError("The simulations directory does not contain any FITS cube files")

    all_measurements_cubes = []
    columns = None

    for cube_path in cube_files:

        print(f"{GREEN}INFO:{RESET} {BOLD}Processing simulation {cube_path.name}{RESET}")
        measurement_sim_cube, col_names = process_simulated_cube(
                cube_path,
                wavelength_range,
                data_extension,
                config_parameters,
                redshift,
                line_restframe,
                real_cube_measured, 
                snr_table, 
                pow)

        print(f"{GREEN}INFO:{RESET} {BOLD}Finished processing simulation {cube_path.name}{RESET}")
        all_measurements_cubes.append(measurement_sim_cube)
        if columns is None:
            columns = col_names
        if columns != col_names:
            raise ValueError("Column mismatch between simulations!")


    all_measurements_cubes = np.stack(all_measurements_cubes)

    print(f"{GREEN}INFO:{RESET} All simulations processed.")

    simulation_results_file = output_dir_path / "simulated_measurements.npy"
    simulation_results_file_npz = output_dir_path / "simulated_measurements.npz"
    print("Saving the array with all measurements from all simulations to", simulation_results_file)

    np.save(simulation_results_file, all_measurements_cubes)
    np.savez(
        simulation_results_file_npz,
        measurements=all_measurements_cubes,
        columns=columns,
        n_sim=len(all_measurements_cubes),
    )
    with open(output_dir_path / "measurements_sim_columns.txt", "w") as f:
        f.write(",".join(columns))

    return all_measurements_cubes