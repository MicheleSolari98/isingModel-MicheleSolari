#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 13:19:36 2026

@author: michele_mac
"""

import sys
from pathlib import Path

from plotter import plot_lattice, plot_magnetization
from storage import load_single_run



def find_latest_single_run(results_directory):
    """Return the most recently saved single-run directory."""
    run_directories = [
        path
        for path in results_directory.glob("single_run_*")
        if path.is_dir()
    ]

    if not run_directories:
        raise FileNotFoundError(
            "No saved single runs were found"
        )

    return max(
        run_directories,
        key=lambda path: path.name,
    )



if __name__ == "__main__":
    results_directory = (
        Path(__file__).parent
        / "results"
    )

    if len(sys.argv) == 1:
        run_directory = find_latest_single_run(
            results_directory
        )

    elif len(sys.argv) == 2:
        run_directory = (
            results_directory
            / sys.argv[1]
        )

    else:
        raise SystemExit(
            "Usage: python plot_single_sim.py "
            "[single_run_directory]"
        )


    results = load_single_run(run_directory)

    sampled_cycles = results["sampled_cycles"]
    magnetization_history = results["magnetization_history"]
    final_lattice = results["final_lattice"]

    if len(magnetization_history) > 0:
        plot_magnetization(
            sampled_cycles,
            magnetization_history,
        )

    plot_lattice(final_lattice)