#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 17:56:51 2026

@author: michele_mac
"""

import numpy as np
from pathlib import Path
import pytest

from ising_simulation import (
    create_lattice,
    magnetization,
    spin_flip_energy,
    run_simulation,
)

from launch_single_sim import (
    load_config,
    validate_single_config,
)

from launch_multi_sim import (
    validate_multi_config,
    create_temperature_grid,
    run_multi_simulation,
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



def test_load_config(tmp_path):
    """A TOML configuration file should be loaded correctly."""
    config_file = tmp_path / "config.toml"

    config_file.write_text(
        """
[physics]
lattice_size = 30
coupling = 1.0
temperature = 2.3
"""
    )

    config = load_config(config_file)

    assert config["physics"]["lattice_size"] == 30
    assert np.isclose(config["physics"]["coupling"], 1.0)
    assert np.isclose(config["physics"]["temperature"], 2.3)


def test_validate_single_config_rejects_invalid_parameters():
    """Invalid simulation parameters should raise a ValueError."""
    invalid_configs = [
        {
            "physics": {
                "lattice_size": 13.7,
                "coupling": 1.0,
                "temperature": 2.3,
            },
            "simulation": {
                "equilibration_cycles": 100,
                "measurement_cycles": 500,
                "seed": 42,
            },
            "output": {
                "record_time_series": True,
                "sample_every": 1,
            },
        },
        {
            "physics": {
                "lattice_size": 30,
                "coupling": 1.0,
                "temperature": -1.0,
            },
            "simulation": {
                "equilibration_cycles": 100,
                "measurement_cycles": 500,
                "seed": 42,
            },
            "output": {
                "record_time_series": True,
                "sample_every": 1,
            },
        },
        {
            "physics": {
                "lattice_size": 30,
                "coupling": 1.0,
                "temperature": 2.3,
            },
            "simulation": {
                "equilibration_cycles": -10,
                "measurement_cycles": 500,
                "seed": 42,
            },
            "output": {
                "record_time_series": True,
                "sample_every": 1,
            },
        },
        {
            "physics": {
                "lattice_size": 30,
                "coupling": 1.0,
                "temperature": 2.3,
            },
            "simulation": {
                "equilibration_cycles": 100,
                "measurement_cycles": 500,
                "seed": 42,
            },
            "output": {
                "record_time_series": True,
                "sample_every": 0,
            },
        },
    ]

    for config in invalid_configs:
        with pytest.raises(ValueError):
            validate_single_config(config)



def test_simulation_sampling_cycles():
    """Time-series samples should be recorded at the requested cycle interval."""
    rng = np.random.default_rng(42)
    lattice = create_lattice(4, rng)

    final_lattice, sampled_cycles, magnetization_history, magnetization_measurements = run_simulation(
        lattice=lattice,
        temperature=2.0,
        coupling=1.0,
        equilibration_cycles=2,
        measurement_cycles=3,
        rng=rng,
        record_time_series=True,
        sample_every=2,
    )

    assert np.array_equal(sampled_cycles, [0, 2, 4])
    assert len(magnetization_history) == len(sampled_cycles)




def test_simulation_is_reproducible():
    """Equal seeds and parameters should produce equal simulations."""
    rng_1 = np.random.default_rng(42)
    rng_2 = np.random.default_rng(42)

    lattice_1 = create_lattice(5, rng_1)
    lattice_2 = create_lattice(5, rng_2)

    final_1, cycles_1, magnetization_1, magnetization_meas1 = run_simulation(
        lattice=lattice_1,
        temperature=2.0,
        coupling=1.0,
        equilibration_cycles=2,
        measurement_cycles=5,
        rng=rng_1,
        record_time_series=True,
        sample_every=1,
    )

    final_2, cycles_2, magnetization_2, magnetization_meas2 = run_simulation(
        lattice=lattice_2,
        temperature=2.0,
        coupling=1.0,
        equilibration_cycles=2,
        measurement_cycles=5,
        rng=rng_2,
        record_time_series=True,
        sample_every=1,
    )

    assert np.array_equal(final_1, final_2)
    assert np.array_equal(cycles_1, cycles_2)
    assert np.array_equal(magnetization_1, magnetization_2)
    
    
    
    
def test_measurement_sampling():
    """Magnetization should be sampled only at the requested measurement interval."""
    rng = np.random.default_rng(42)
    lattice = create_lattice(4, rng)

    (
        final_lattice,
        sampled_cycles,
        magnetization_history,
        magnetization_measurements,
    ) = run_simulation(
        lattice=lattice,
        temperature=2.0,
        coupling=1.0,
        equilibration_cycles=2,
        measurement_cycles=10,
        rng=rng,
        record_time_series=False,
        sample_every=1,
        measurement_every=2,
    )

    assert len(magnetization_measurements) == 5
    
    
    
def test_temperature_grid_includes_critical_temperature():
    """The temperature grid should include Tc for J = 1 when requested."""
    temperatures = create_temperature_grid(
        start=1.0,
        stop=3.5,
        points=5,
        coupling=1.0,
        include_critical_temperature=True,
    )

    assert np.any(
        np.isclose(
            temperatures,
            2.269,
            rtol=0,
            atol=1e-3,
        )
    )
    
    
    
    
def test_multi_run_output_shape():
    """The multi-run launcher should return the expected result shape."""
    config = {
        "physics": {
            "lattice_size": 4,
            "coupling": 1.0,
        },
        "simulation": {
            "equilibration_cycles": 1,
            "measurement_cycles": 4,
            "measurement_every": 2,
            "seed": 42,
            "repetitions": 2,
        },
        "temperature_sweep": {
            "start": 1.0,
            "stop": 3.0,
            "points": 3,
            "include_critical_temperature": True,
        },
        "output": {
            "record_time_series": False,
            "sample_every": 1,
        },
    }

    temperatures, magnetization_measurements = (
        run_multi_simulation(config)
    )

    assert len(temperatures) == 4
    assert magnetization_measurements.shape == (4, 2, 2)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    