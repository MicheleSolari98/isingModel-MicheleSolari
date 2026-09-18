#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 13:11:18 2026

@author: michele_mac
"""

from datetime import datetime
from pathlib import Path
import shutil

import numpy as np


def save_single_run(
    results_directory,
    config_path,
    sampled_cycles,
    magnetization_history,
    final_lattice,
):
    """Save the raw results and configuration of a single simulation."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    run_directory = (
        Path(results_directory)
        / f"single_run_{timestamp}"
    )

    run_directory.mkdir(
        parents=True,
        exist_ok=False,
    )

    np.savez_compressed(
        run_directory / "data.npz",
        sampled_cycles=sampled_cycles,
        magnetization_history=magnetization_history,
        final_lattice=final_lattice,
    )

    shutil.copy2(
        config_path,
        run_directory / "parameters.toml",
    )

    return run_directory


def load_single_run(run_directory):
    """Load the raw results of a saved single simulation."""
    data_path = Path(run_directory) / "data.npz"

    with np.load(data_path) as data:
        results = {
            key: data[key]
            for key in data.files
        }

    return results