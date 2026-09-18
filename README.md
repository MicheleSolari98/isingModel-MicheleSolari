# 2D Ising Model Monte Carlo Simulation

Monte Carlo simulation of the two-dimensional Ising model on a square lattice using the Metropolis algorithm.

The project supports both single-temperature simulations and temperature sweeps, together with statistical analysis, comparison with exact analytical results, data storage, and reproducible plotting.

## Features

- 2D square Ising lattice with periodic boundary conditions
- Metropolis Monte Carlo dynamics
- Single-temperature simulations
- Temperature sweeps with independent repetitions
- Magnetization measurements and statistical uncertainties
- Energy calculated from final configurations
- Comparison with exact thermodynamic-limit results
- Saving and reloading simulation results
- Reproducible random-number generation through explicit seeds
- Automated tests with `pytest`

## Model

The Hamiltonian of the zero-field Ising model is

$$
H = -J \sum_{\langle i,j \rangle} s_i s_j,
$$

where

- $s_i = \pm 1$ is the spin at lattice site $i$,
- $J$ is the nearest-neighbor coupling,
- the sum runs over nearest-neighbor pairs.

The simulations use units where

$$
k_B = 1.
$$

For the ferromagnetic model, the exact critical temperature is

$$
T_c =
\frac{2J}
{\ln(1+\sqrt{2})}.
$$

For $J=1$,

$$
T_c \approx 2.269.
$$

## Project structure

```text
.
├── configs/
│   ├── single_run_parameters.toml
│   └── multi_run_parameters.toml
├── results/
├── analysis.py
├── analyze_multi_sim.py
├── ising_simulation.py
├── launch_multi_sim.py
├── launch_single_sim.py
├── plot_single_sim.py
├── plotter.py
├── storage.py
└── test_simulation.py