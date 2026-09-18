#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 22:36:26 2026

@author: michele_mac
"""

import matplotlib.pyplot as plt


def plot_magnetization(sampled_cycles, magnetization_history):
    """Plot magnetization as a function of simulation cycles."""
    plt.plot(sampled_cycles, magnetization_history)
    plt.xlabel("Simulation cycle")
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
    
    
    
    
def plot_magnetization_vs_temperature(
    temperatures,
    measured_magnetization,
    magnetization_error,
    theoretical_temperatures,
    theoretical_magnetization,
):
    """Plot measured and theoretical magnetization versus temperature."""
    plt.errorbar(
        temperatures,
        measured_magnetization,
        yerr=magnetization_error,
        fmt="o",
        capsize=3,
        label="Monte Carlo <|m|>",
    )

    plt.plot(
        theoretical_temperatures,
        theoretical_magnetization,
        label="Exact solution",
    )

    plt.xlabel("Temperature")
    plt.ylabel("Absolute magnetization per site")
    plt.title("Magnetization versus temperature")
    plt.legend()
    plt.show()
    
    
    
def plot_energy_vs_temperature(
    temperatures,
    measured_energy,
    theoretical_temperatures,
    theoretical_energy,
):
    """Plot measured and theoretical energy versus temperature."""
    plt.plot(
        temperatures,
        measured_energy,
        "o",
        label="Monte Carlo final-state energy",
    )

    plt.plot(
        theoretical_temperatures,
        theoretical_energy,
        label="Exact energy",
    )

    plt.xlabel("Temperature")
    plt.ylabel("Energy per site")
    plt.title("Energy versus temperature")
    plt.legend()
    plt.show()
    
    