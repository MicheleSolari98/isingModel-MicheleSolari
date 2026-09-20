#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 14:01:43 2026

@author: michele_mac
"""

import sys
import tomllib
from pathlib import Path

import numpy as np

from analysis import (
    magnetization_statistics,
    onsager_magnetization,
    mean_final_energy,
    onsager_energy,
    mean_connected_correlation,
)

from plotter import (
    plot_magnetization_vs_temperature, 
    plot_energy_vs_temperature,
    plot_connected_correlation_vs_temperature
)

from storage import load_multi_run



def find_latest_multi_run(results_directory):
    """Return the most recently saved multi-run directory."""
    run_directories = [
        path
        for path in results_directory.glob("multi_run_*")
        if path.is_dir()
    ]

    if not run_directories:
        raise FileNotFoundError(
            "No saved multi runs were found"
        )

    return max(
        run_directories,
        key=lambda path: path.name,
    )


def analyze_and_plot_multi_run(results, config):
    """Analyze and plot the magnetization of a multi run."""
    temperatures = results["temperatures"]
    magnetization_measurements = results["magnetization"]
    
    final_lattices = results["final_lattices"]
    coupling = config["physics"]["coupling"]
    # Distance used for the connected-correlation plot
    correlation_distance_divisor = 10

    lattice_size = final_lattices.shape[-1]

    correlation_distance = max(
    1,
    lattice_size // correlation_distance_divisor,
    )

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
        temperatures.max()
        - temperatures.min()
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
    
    #Magnetization
    theoretical_magnetization = onsager_magnetization(
        theoretical_temperatures,
        coupling=coupling,
    )

    plot_magnetization_vs_temperature(
        temperatures,
        mean_magnetization,
        magnetization_error,
        theoretical_temperatures,
        theoretical_magnetization,
    )
    
    
    # Energy
    mean_energy = mean_final_energy(
        final_lattices,
        coupling=coupling,
    )

    theoretical_energy = onsager_energy(
        theoretical_temperatures,
        coupling=coupling,
    )

    plot_energy_vs_temperature(
        temperatures,
        mean_energy,
        theoretical_temperatures,
        theoretical_energy,
    )
    
    # Connected spin correlation
    mean_correlation = mean_connected_correlation(
        final_lattices,
        distance=correlation_distance,
    )

    plot_connected_correlation_vs_temperature(
        temperatures,
        mean_correlation,
        correlation_distance,
    )







if __name__ == "__main__":
    results_directory = (
        Path(__file__).parent
        / "results"
    )

    if len(sys.argv) == 1:
        run_directory = find_latest_multi_run(
            results_directory
        )

    elif len(sys.argv) == 2:
        run_directory = (
            results_directory
            / sys.argv[1]
        )

    else:
        raise SystemExit(
            "Usage: python analyze_multi_sim.py "
            "[multi_run_directory]"
        )
        
        
    results = load_multi_run(
        run_directory
    )

    config_path = (
        run_directory
        / "parameters.toml"
    )

    with open(config_path, "rb") as file:
        config = tomllib.load(file)

    print(
        "Analyzing results from:",
        run_directory.name,
    )

    analyze_and_plot_multi_run(
        results,
        config,
    )