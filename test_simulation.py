#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 17:56:51 2026

@author: michele_mac
"""

import numpy as np

from ising_simulation import (
    create_lattice,
    magnetization,
    spin_flip_energy,
    run_simulation,
)


def test_create_lattice_contains_valid_spins():
    """The lattice should contain only integer spins -1 and +1."""
    rng = np.random.default_rng(42)

    lattice = create_lattice(10, rng)

    assert np.issubdtype(lattice.dtype, np.integer)
    assert set(np.unique(lattice)).issubset({-1, 1})


def test_magnetization_known_configuration():
    """Magnetization should be correct for a configuration with known result."""
    lattice = np.array([
        [1, 1],
        [-1, -1],
    ])

    assert np.isclose(magnetization(lattice), 0.0)


def test_flip_energy_aligned_spin():
    """Flipping a spin aligned with four equal neighbours should cost 8J."""
    lattice = np.ones((3, 3), dtype=int)

    delta_energy = spin_flip_energy(
        lattice,
        row=1,
        col=1,
        coupling=1.0,
    )

    assert np.isclose(delta_energy, 8.0)


def test_simulation_is_reproducible():
    """Equal seeds and parameters should produce equal simulations."""
    rng_1 = np.random.default_rng(42)
    rng_2 = np.random.default_rng(42)

    lattice_1 = create_lattice(5, rng_1)
    lattice_2 = create_lattice(5, rng_2)

    final_1, magnetization_1 = run_simulation(
        lattice_1, 2.0, 1.0, 10, rng_1
    )
    final_2, magnetization_2 = run_simulation(
        lattice_2, 2.0, 1.0, 10, rng_2
    )

    assert np.array_equal(final_1, final_2)
    assert np.array_equal(magnetization_1, magnetization_2)