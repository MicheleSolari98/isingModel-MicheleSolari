#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 00:19:32 2026

@author: michele_mac
"""

import tomllib
from pathlib import Path

import numpy as np

from ising_simulation import create_lattice, run_simulation
from plotter import plot_lattice, plot_magnetization




def load_config(config_path):
    """Load simulation parameters from a TOML configuration file."""
    with open(config_path, "rb") as file:
        return tomllib.load(file)



def validate_single_config(config):
    """Validate the parameters loaded from the configuration file."""
    try:
        lattice_size = config["physics"]["lattice_size"]
        coupling = config["physics"]["coupling"]
        temperature = config["physics"]["temperature"]

        equilibration_cycles = config["simulation"]["equilibration_cycles"]
        measurement_cycles = config["simulation"]["measurement_cycles"]
        seed = config["simulation"]["seed"]

        record_time_series = config["output"]["record_time_series"]
        sample_every = config["output"]["sample_every"]

    except (KeyError, TypeError) as error:
        raise ValueError("Missing or invalid configuration structure") from error

    if type(lattice_size) is not int or lattice_size <= 0:
        raise ValueError("lattice_size must be a positive integer")

    if (
        not isinstance(coupling, (int, float))
        or isinstance(coupling, bool)
    ):
        raise ValueError("coupling must be a number")

    if (
        not isinstance(temperature, (int, float))
        or isinstance(temperature, bool)
        or temperature <= 0
    ):
        raise ValueError("temperature must be a positive number")

    if (
        type(equilibration_cycles) is not int
        or equilibration_cycles < 0
    ):
        raise ValueError(
            "equilibration_cycles must be a non-negative integer"
        )

    if (
        type(measurement_cycles) is not int
        or measurement_cycles < 0
    ):
        raise ValueError(
            "measurement_cycles must be a non-negative integer"
        )

    if type(seed) is not int or seed < 0:
        raise ValueError("seed must be a non-negative integer")

    if type(record_time_series) is not bool:
        raise ValueError("record_time_series must be a boolean")

    if type(sample_every) is not int or sample_every <= 0:
        raise ValueError("sample_every must be a positive integer")











if __name__ == "__main__":
    config_path = (
        Path(__file__).parent
        / "configs"
        / "single_run_parameters.toml"
    )

    config = load_config(config_path)
    validate_single_config(config)

    # Physical parameters
    lattice_size = config["physics"]["lattice_size"]
    coupling = config["physics"]["coupling"]
    temperature = config["physics"]["temperature"]

    # Simulation parameters
    equilibration_cycles = config["simulation"]["equilibration_cycles"]
    measurement_cycles = config["simulation"]["measurement_cycles"]
    seed = config["simulation"]["seed"]

    # Output parameters
    record_time_series = config["output"]["record_time_series"]
    sample_every = config["output"]["sample_every"]

    # Random number generator
    rng = np.random.default_rng(seed)

    # Initial random spin configuration
    lattice = create_lattice(lattice_size, rng)

    # Run simulation
    (
        final_lattice,
        sampled_cycles,
        magnetization_history,
        magnetization_measurements,
    ) = run_simulation(
        lattice=lattice,
        temperature=temperature,
        coupling=coupling,
        equilibration_cycles=equilibration_cycles,
        measurement_cycles=measurement_cycles,
        rng=rng,
        record_time_series=record_time_series,
        sample_every=sample_every,
    )

    # Simple visual checks
    if record_time_series:
        plot_magnetization(
            sampled_cycles,
            magnetization_history,
        )

    plot_lattice(final_lattice)








