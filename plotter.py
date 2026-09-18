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
    theoretical_temperatures,
    theoretical_magnetization,
):
    """Plot measured and theoretical magnetization versus temperature."""
    plt.plot(
        temperatures,
        measured_magnetization,
        "o",
        label="Monte Carlo",
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