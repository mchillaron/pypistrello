#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

"""Read a YAML file and extract information from it"""

import yaml

def load_yaml_file(config_path):
    """
    Load YAML configuration file into a Python dictionary.
    """
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    if config is None:
        raise ValueError("Configuration file is empty.")

    return config

    
def validate_region_config(config):
    """
    Validate mutual exclusivity and required fields.
    """

    # Line region
    if config["line"].get("window_line") is not None and config["line"].get("fit_region") is not None:
        raise ValueError("Both window_line and fit_region are defined. Choose only one.")

    if config["line"].get("window_line") is None and config["line"].get("fit_region") is None:
        raise ValueError("Neither window_line nor fit_region is defined.")

    # Continuum region
    if config["continuum"].get("window_cont") is not None and config["continuum"].get("continuum_region") is not None:
        raise ValueError("Both window_cont and continuum_region are defined. Choose only one.")

    if config["continuum"].get("window_cont") is None and config["continuum"].get("continuum_region") is None:
        raise ValueError("Neither window_cont nor continuum_region is defined.")
