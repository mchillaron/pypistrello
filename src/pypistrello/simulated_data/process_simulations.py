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

def process_simulations(
        simulation_dir,
        output_dir_path,
        wavelength_range,
        data_extension,
        config_parameters,
        redshift,
        line_restframe):

    cube_files = sorted(
        f for f in simulation_dir.iterdir()
        if f.suffix.lower() in [".fits", ".fit"]
    )
    print(f"{GREEN}INFO:{RESET} A total number of {len(cube_files)} simulation cubes will be processed")

    all_measurements_areatrapz_cubes = []

    for cube_path in cube_files:

        print(f"{GREEN}INFO:{RESET} Processing simulation {cube_path.name}")

        measurement_areatrapz_sim_cube = process_simulated_cube(
            cube_path,
            wavelength_range,
            data_extension,
            config_parameters,
            redshift,
            line_restframe
        )

        print(f"{GREEN}INFO:{RESET} Finished processing simulation {cube_path.name}")
        all_measurements_areatrapz_cubes.append(measurement_areatrapz_sim_cube)

    all_measurements_areatrapz_cubes = np.stack(all_measurements_areatrapz_cubes)
    #print(all_measurements_areatrapz_cubes)

    print(f"{GREEN}INFO:{RESET} All simulations processed.")

    simulation_results_file = output_dir_path / "simulated_measurements.npy"
    simulation_results_file_npz = output_dir_path / "simulated_measurements.npz"
    print("Saving the array with all measurements from all simulations to", simulation_results_file)

    np.save(simulation_results_file, all_measurements_areatrapz_cubes)
    measurement_names = ["x", "y", "flux", "noise_cont", "snr_cont"]
    with open("measuremens_sim_columns.txt", "w") as f:
        f.write(",".join(measurement_names))

    np.savez(
        simulation_results_file_npz, 
        measurements=all_measurements_areatrapz_cubes,
        columns=measurement_names
    )

    return all_measurements_areatrapz_cubes