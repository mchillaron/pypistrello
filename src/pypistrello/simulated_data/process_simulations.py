#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from astropy.table import Table

import numpy as np

from .process_simulated_cube import process_simulated_cube

GREEN   = "\033[92m"
RESET   = "\033[0m"
BOLD = "\033[1m"



def array_to_table(array, columns):
    """
    Reconstruct an Astropy Table from a numpy array and a list of column names.

    Parameters
    ----------
    array : ndarray
        Shape (N_spaxels, N_columns)

    columns : list
        Column names

    Returns
    -------
    table : astropy.table.Table
    """

    COLUMN_DTYPES = {
        "x": int,
        "y": int,
        "bin_id": int,
        "n_pix": int,
        "area_trapz": float,
        "cont_noise": float,
        "snr_trapz": float,
        "velocity": float,
        "offsets": float,
    }
    
    table = Table()

    #for i, col in enumerate(columns):
    #    #table[col] = array[:, i]
    #    table[col] = array[:, i].copy()

    for i, col in enumerate(columns):

        dtype = COLUMN_DTYPES.get(col)

        if dtype is not None:
            table[col] = array[:, i].astype(dtype)
        else:
            table[col] = array[:, i]

    return table



def process_simulations(
        simulation_dir,
        output_dir_path,
        wavelength_range,
        data_extension,
        config_parameters,
        redshift,
        line_restframe,
        sim_results_file,
        real_cube_measured=False, 
        trapz_npz_filename=None,
        pow=None,
        pow_valid_mask=None,
        ):

    cube_files = sorted(
        f for f in simulation_dir.iterdir()
        if f.suffix.lower() in [".fits", ".fit"]
    )
    print(f"{GREEN}INFO:{RESET} A total number of {len(cube_files)} simulation cubes will be processed")
    if len(cube_files) == 0:
        ValueError("The simulations directory does not contain any FITS cube files")

    if real_cube_measured and trapz_npz_filename is not None:
        print(f"{GREEN}INFO:{RESET} Loading trapezoidal measurements from {trapz_npz_filename}")
        data_trapz_npz = np.load(trapz_npz_filename, allow_pickle=True)

        measurements_trapz = data_trapz_npz["measurements"]
        print(f"{GREEN}INFO:{RESET} Shape of trapezoidal measurements: {measurements_trapz.shape}")
        columns_trapz = data_trapz_npz["columns"].tolist()

        print()
        print(f"INFO: Loaded {measurements_trapz.shape[0]} simulations")
        
        if measurements_trapz.shape[0] != len(cube_files):
            raise ValueError(
                "The number of saved simulations does not match "
                "the number of FITS cubes."
            )

    else:
        trapz_table = None

    all_measurements_cubes = []
    columns = None

    #for cube_path in cube_files:
    for i, cube_path in enumerate(cube_files):

        if real_cube_measured:
            trapz_table = array_to_table(measurements_trapz[i], columns_trapz)

        print(f"{GREEN}INFO:{RESET} {BOLD}Processing simulation {cube_path.name}{RESET}")
        measurement_sim_cube, col_names = process_simulated_cube(
                cube_path,
                wavelength_range,
                data_extension,
                config_parameters,
                redshift,
                line_restframe,
                real_cube_measured, 
                trapz_table,
                pow,
                pow_valid_mask)

        print(f"{GREEN}INFO:{RESET} {BOLD}Finished processing simulation {cube_path.name}{RESET}")
        print()
        
        all_measurements_cubes.append(measurement_sim_cube)
        if columns is None:
            columns = col_names
        if columns != col_names:
            raise ValueError("Column mismatch between simulations!")


    all_measurements_cubes = np.stack(all_measurements_cubes)

    print(f"{GREEN}INFO:{RESET} All simulations processed.")

    
    print("Saving the array with all measurements from all simulations to", sim_results_file)

    #np.save(simulation_results_file, all_measurements_cubes)
    np.savez(
        sim_results_file,
        measurements=all_measurements_cubes,
        columns=columns,
        n_sim=len(all_measurements_cubes),
    )
    with open(output_dir_path / "measurements_sim_columns.txt", "w") as f:
        f.write(",".join(columns))

    return all_measurements_cubes