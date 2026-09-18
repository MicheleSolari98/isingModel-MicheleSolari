#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 02:06:27 2026

@author: michele_mac
"""

import numpy as np


def mean_absolute_magnetization(magnetization_measurements):
    """Calculate mean absolute magnetization for each temperature."""
    return np.mean(
        np.abs(magnetization_measurements),
        axis=(1, 2),
    )


def onsager_magnetization(temperatures, coupling):
    """Calculate the exact spontaneous magnetization of the 2D Ising model."""
    if coupling <= 0:
        raise ValueError(
            "Onsager magnetization comparison requires positive coupling"
        )

    critical_temperature = (
        2 * coupling
        / np.log(1 + np.sqrt(2))
    )

    theoretical_magnetization = np.zeros_like(
        temperatures,
        dtype=float,
    )

    below_critical = temperatures < critical_temperature

    theoretical_magnetization[below_critical] = (
        1
        - np.sinh(
            2 * coupling / temperatures[below_critical]
        ) ** (-4)
    ) ** (1 / 8)

    return theoretical_magnetization