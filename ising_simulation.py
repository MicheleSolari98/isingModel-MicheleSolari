#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 14:46:13 2023

@author: michele_mac
"""
import math
import matplotlib.pyplot as plt
import numpy as np
import tomllib
from pathlib import Path


def load_config(config_path):
    """Load simulation parameters from a TOML configuration file."""
    with open(config_path, "rb") as file:
        return tomllib.load(file)


def validate_config(config):
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




def create_lattice(size, rng):
    """Create a square lattice with randomly oriented spins (+1 or -1)."""
    return rng.choice([-1, 1], size=(size, size))


def magnetization(lattice):
    """Return the mean magnetization per lattice site."""
    return np.mean(lattice)


def spin_flip_energy(lattice, row, col, coupling):
    """Calculate the energy change caused by flipping one spin."""
    size = lattice.shape[0]
    spin = lattice[row, col]

    neighbors = (
        lattice[(row - 1) % size, col]
        + lattice[(row + 1) % size, col]
        + lattice[row, (col - 1) % size]
        + lattice[row, (col + 1) % size]
    )

    return 2 * coupling * spin * neighbors


def metropolis_step(lattice, temperature, coupling, rng):
    """Attempt one spin flip according to the Metropolis algorithm."""
    size = lattice.shape[0]

    row = rng.integers(0, size)
    col = rng.integers(0, size)

    delta_energy = spin_flip_energy(
        lattice,
        row,
        col,
        coupling,
    )

    if delta_energy <= 0:
        lattice[row, col] *= -1

    elif rng.random() < math.exp(-delta_energy / temperature):
        lattice[row, col] *= -1

def run_cycle(lattice, temperature, coupling, rng):
    """Run one simulation cycle, corresponding to L^2 flip attempts."""
    for i in range(lattice.size):
        metropolis_step(
            lattice,
            temperature,
            coupling,
            rng,
        )



def run_simulation(
    lattice,
    temperature,
    coupling,
    equilibration_cycles,
    measurement_cycles,
    rng,
    record_time_series,
    sample_every,
):
    """
    Run the Ising simulation.

    The system is first equilibrated and then evolved during the
    measurement phase. If requested, magnetization is recorded
    periodically throughout the simulation.
    """
    sampled_cycles = []
    magnetization_history = []

    if record_time_series:
        sampled_cycles.append(0)
        magnetization_history.append(magnetization(lattice))

    # Equilibration phase
    for i in range(1, equilibration_cycles + 1):
        run_cycle(
            lattice,
            temperature,
            coupling,
            rng,
        )

        if record_time_series and i % sample_every == 0:
            sampled_cycles.append(i)
            magnetization_history.append(magnetization(lattice))

    # Measurement phase
    for i in range(1, measurement_cycles + 1):
        run_cycle(
            lattice,
            temperature,
            coupling,
            rng,
        )

        current_cycle = equilibration_cycles + i

        if record_time_series and current_cycle % sample_every == 0:
            sampled_cycles.append(current_cycle)
            magnetization_history.append(magnetization(lattice))

    return (
        lattice,
        np.array(sampled_cycles),
        np.array(magnetization_history),
    )


def plot_magnetization(sampled_cycles, magnetization_history):
    """Plot magnetization as a function of Monte Carlo sweeps."""
    plt.plot(sampled_cycles, magnetization_history)
    plt.xlabel("Monte Carlo sweep")
    plt.ylabel("Magnetization per site")
    plt.title("Magnetization evolution")
    plt.show()


def plot_lattice(lattice):
    """Plot the spin configuration of the lattice."""
    plt.imshow(
        lattice,
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        interpolation="nearest",
    )

    plt.title("Final spin configuration")
    plt.axis("off")
    plt.show()


if __name__ == "__main__":
    config_path = (
        Path(__file__).parent
        / "configs"
        / "single_run_parameters.toml"
    )

    config = load_config(config_path)
    validate_config(config)

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
