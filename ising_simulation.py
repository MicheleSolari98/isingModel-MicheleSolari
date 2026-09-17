#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 14:46:13 2023

@author: michele_mac
"""
import math
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
    measurement_every=None,
):
    """
    Run the Ising simulation.

    The system is first equilibrated and then evolved during the
    measurement phase. If requested, magnetization is recorded
    periodically throughout the simulation.
    """
    sampled_cycles = []
    magnetization_history = []
    magnetization_measurements = []

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
        if (
            measurement_every is not None
            and i % measurement_every == 0
        ):
            magnetization_measurements.append(magnetization(lattice))

    return (
        lattice,
        np.array(sampled_cycles),
        np.array(magnetization_history),
        np.array(magnetization_measurements),
    )





