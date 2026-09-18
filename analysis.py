#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 02:06:27 2026

@author: michele_mac
"""

import numpy as np


def autocorrelation_factor(values):
    """Estimate the autocorrelation factor of a measurement series."""
    centered_values = values - np.mean(values)

    variance_sum = np.dot(
        centered_values,
        centered_values,
    )

    if variance_sum == 0:
        return 1.0

    g = 1.0

    for lag in range(1, len(values)):
        autocorrelation = (
            np.dot(
                centered_values[:-lag],
                centered_values[lag:],
            )
            / variance_sum
        )

        if autocorrelation <= 0:
            break

        g += 2 * autocorrelation

    return g




def magnetization_statistics(magnetization_measurements):
    """
    Calculate mean absolute magnetization and its statistical error.

    Measurements from each independent run are first analyzed separately.
    Autocorrelation is used to estimate the effective number of independent
    measurements. The final uncertainty also includes additional variance
    observed between independent runs (different seeds).
    """
    absolute_measurements = np.abs(
        magnetization_measurements
    )

    repetitions = absolute_measurements.shape[1]
    number_of_measurements = absolute_measurements.shape[2]

    if number_of_measurements < 2:
        raise ValueError(
            "At least two measurements per run are required"
        )

    mean_per_seed = np.mean(
        absolute_measurements,
        axis=2,
    )

    error_per_seed = np.empty_like(
        mean_per_seed,
        dtype=float,
    )

    for i in range(absolute_measurements.shape[0]):
        for j in range(repetitions):
            measurements = absolute_measurements[i, j]

            g = autocorrelation_factor(measurements)

            effective_measurements = (
                number_of_measurements / g
            )

            standard_deviation = np.std(
                measurements,
                ddof=1,
            )

            error_per_seed[i, j] = (
                standard_deviation
                / np.sqrt(effective_measurements)
            )

    mean_magnetization = np.mean(
        mean_per_seed,
        axis=1,
    )

    magnetization_error = np.empty(
        absolute_measurements.shape[0]
    )

    for i in range(absolute_measurements.shape[0]):
        if repetitions == 1:
            magnetization_error[i] = error_per_seed[i, 0]

        else:
            between_seed_variance = np.var(
                mean_per_seed[i],
                ddof=1,
            )

            mean_within_variance = np.mean(
                error_per_seed[i] ** 2
            )

            extra_between_variance = max(
                0.0,
                between_seed_variance
                - mean_within_variance,
            )

            within_variance_of_mean = (
                np.sum(error_per_seed[i] ** 2)
                / repetitions**2
            )

            magnetization_error[i] = np.sqrt(
                extra_between_variance / repetitions
                + within_variance_of_mean
            )

    return mean_magnetization, magnetization_error





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