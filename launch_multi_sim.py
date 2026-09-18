#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 01:27:12 2026

@author: michele_mac
"""

import tomllib
from pathlib import Path

import numpy as np

from ising_simulation import create_lattice, run_simulation

from analysis import (
    magnetization_statistics,
    onsager_magnetization,
)

from plotter import plot_magnetization_vs_temperature



def load_config(config_path):
    """Load simulation parameters from a TOML configuration file."""
    with open(config_path, "rb") as file:
        return tomllib.load(file)
    
    
    
def validate_multi_config(config):
    """Validate the multi-run configuration parameters."""
    try:
        lattice_size = config["physics"]["lattice_size"]
        coupling = config["physics"]["coupling"]

        equilibration_cycles = config["simulation"]["equilibration_cycles"]
        measurement_cycles = config["simulation"]["measurement_cycles"]
        measurement_every = config["simulation"]["measurement_every"]
        seed = config["simulation"]["seed"]
        repetitions = config["simulation"]["repetitions"]

        temperature_start = config["temperature_sweep"]["start"]
        temperature_stop = config["temperature_sweep"]["stop"]
        temperature_points = config["temperature_sweep"]["points"]
        include_critical_temperature = (
            config["temperature_sweep"]["include_critical_temperature"]
        )

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
        type(equilibration_cycles) is not int
        or equilibration_cycles < 0
    ):
        raise ValueError(
            "equilibration_cycles must be a non-negative integer"
        )

    if (
        type(measurement_cycles) is not int
        or measurement_cycles <= 0
    ):
        raise ValueError(
            "measurement_cycles must be a positive integer"
        )

    if type(measurement_every) is not int or measurement_every <= 0:
        raise ValueError("measurement_every must be a positive integer")

    if measurement_every > measurement_cycles:
        raise ValueError(
            "measurement_every cannot exceed measurement_cycles"
        )
    
    if measurement_cycles // measurement_every < 2:
        raise ValueError(
            "At least two measurements per run are required"
        )
    
    
    if type(seed) is not int or seed < 0:
        raise ValueError("seed must be a non-negative integer")

    if type(repetitions) is not int or repetitions <= 0:
        raise ValueError("repetitions must be a positive integer")

    if (
        not isinstance(temperature_start, (int, float))
        or isinstance(temperature_start, bool)
        or temperature_start <= 0
    ):
        raise ValueError("temperature sweep start must be positive")

    if (
        not isinstance(temperature_stop, (int, float))
        or isinstance(temperature_stop, bool)
        or temperature_stop <= temperature_start
    ):
        raise ValueError(
            "temperature sweep stop must be greater than start"
        )

    if type(temperature_points) is not int or temperature_points < 2:
        raise ValueError("temperature sweep points must be at least 2")

    if type(include_critical_temperature) is not bool:
        raise ValueError(
            "include_critical_temperature must be a boolean"
        )

    if include_critical_temperature and coupling == 0:
        raise ValueError(
            "critical temperature requires non-zero coupling"
        )

    if type(record_time_series) is not bool:
        raise ValueError("record_time_series must be a boolean")

    if type(sample_every) is not int or sample_every <= 0:
        raise ValueError("sample_every must be a positive integer")
        
        
        
        
def create_temperature_grid(
    start,
    stop,
    points,
    coupling,
    include_critical_temperature,
):
    """Create the temperature grid for the simulation sweep."""
    temperatures = np.linspace(
        start,
        stop,
        points,
    )

    if include_critical_temperature:
        critical_temperature = (
            2 * abs(coupling)
            / np.log(1 + np.sqrt(2))
        )

        temperatures = np.append(
            temperatures,
            critical_temperature,
        )

        temperatures = np.unique(
            np.sort(temperatures)
        )

    return temperatures





def run_multi_simulation(config):
    """Run simulations over all temperatures and repetitions."""
    validate_multi_config(config)

    # Physical parameters
    lattice_size = config["physics"]["lattice_size"]
    coupling = config["physics"]["coupling"]

    # Simulation parameters
    equilibration_cycles = config["simulation"]["equilibration_cycles"]
    measurement_cycles = config["simulation"]["measurement_cycles"]
    measurement_every = config["simulation"]["measurement_every"]
    seed = config["simulation"]["seed"]
    repetitions = config["simulation"]["repetitions"]

    # Temperature sweep parameters
    temperature_start = config["temperature_sweep"]["start"]
    temperature_stop = config["temperature_sweep"]["stop"]
    temperature_points = config["temperature_sweep"]["points"]
    include_critical_temperature = (
        config["temperature_sweep"]["include_critical_temperature"]
    )

    # Output parameters
    record_time_series = config["output"]["record_time_series"]
    sample_every = config["output"]["sample_every"]

    temperatures = create_temperature_grid(
        start=temperature_start,
        stop=temperature_stop,
        points=temperature_points,
        coupling=coupling,
        include_critical_temperature=include_critical_temperature,
    )

    number_of_measurements = (
        measurement_cycles // measurement_every
    )

    magnetization_measurements = np.empty(
        (
            len(temperatures),
            repetitions,
            number_of_measurements,
        )
    )

    for i in range(len(temperatures)):
        temperature = temperatures[i]

        for j in range(repetitions):
            current_seed = (
                seed
                + i * repetitions
                + j
            )

            rng = np.random.default_rng(current_seed)

            lattice = create_lattice(
                lattice_size,
                rng,
            )

            (
                final_lattice,
                sampled_cycles,
                magnetization_history,
                current_measurements,
            ) = run_simulation(
                lattice=lattice,
                temperature=temperature,
                coupling=coupling,
                equilibration_cycles=equilibration_cycles,
                measurement_cycles=measurement_cycles,
                rng=rng,
                record_time_series=record_time_series,
                sample_every=sample_every,
                measurement_every=measurement_every,
            )

            magnetization_measurements[i, j] = (
                current_measurements
            )

    return temperatures, magnetization_measurements







if __name__ == "__main__":
    config_path = (
        Path(__file__).parent
        / "configs"
        / "multi_run_parameters.toml"
    )

    config = load_config(config_path)

    temperatures, magnetization_measurements = (
        run_multi_simulation(config)
    )

    # Analyze magnetization measurements
    (
    mean_magnetization,
    magnetization_error,
) = magnetization_statistics(
    magnetization_measurements
)
    # Parameters for the theoretical curve
    plot_margin_fraction = 0.05
    theoretical_points = 500

    plot_margin = plot_margin_fraction * (
        temperatures.max() - temperatures.min()
    )

    plot_temperature_min = max(
        temperatures.min() - plot_margin,
        1e-6,
    )

    plot_temperature_max = (
        temperatures.max() + plot_margin
    )

    theoretical_temperatures = np.linspace(
        plot_temperature_min,
        plot_temperature_max,
        theoretical_points,
    )

    theoretical_magnetization = onsager_magnetization(
        theoretical_temperatures,
        coupling=config["physics"]["coupling"],
    )

    # Plot Monte Carlo results and exact solution
    plot_magnetization_vs_temperature(
        temperatures,
        mean_magnetization,
        magnetization_error,
        theoretical_temperatures,
        theoretical_magnetization,
    )









