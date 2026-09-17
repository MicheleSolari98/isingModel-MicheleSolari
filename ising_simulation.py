#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 14:46:13 2023

@author: michele_mac
"""
import math

import matplotlib.pyplot as plt
import numpy as np


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


def run_simulation(
    lattice,
    temperature,
    coupling,
    n_sweeps,
    rng,
):
    """
    Run the Metropolis simulation.

    One Monte Carlo sweep corresponds to L^2 attempted spin flips,
    where L is the lattice size.

    Returns the final lattice and the magnetization measured after
    each sweep.
    """
    size = lattice.shape[0]
    attempts_per_sweep = size**2

    magnetization_history = [magnetization(lattice)]

    for _ in range(n_sweeps):
        for _ in range(attempts_per_sweep):
            metropolis_step(
                lattice,
                temperature,
                coupling,
                rng,
            )

        magnetization_history.append(magnetization(lattice))

    return lattice, np.array(magnetization_history)


def plot_magnetization(magnetization_history):
    """Plot magnetization as a function of Monte Carlo sweeps."""
    sweeps = np.arange(len(magnetization_history))

    plt.plot(sweeps, magnetization_history)
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
    # Simulation parameters
    L = 30
    T = 2.3
    J = 1.0
    N_SWEEPS = 1000
    SEED = 42

    # Random number generator
    rng = np.random.default_rng(SEED)

    # Initial random spin configuration
    lattice = create_lattice(L, rng)

    # Run simulation
    final_lattice, magnetization_history = run_simulation(
        lattice=lattice,
        temperature=T,
        coupling=J,
        n_sweeps=N_SWEEPS,
        rng=rng,
    )

    # Simple visual checks
    plot_magnetization(magnetization_history)
    plot_lattice(final_lattice)

