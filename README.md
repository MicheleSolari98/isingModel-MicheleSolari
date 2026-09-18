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

## Physical model

The zero-field two-dimensional Ising model is defined by the Hamiltonian

$$
H = -J \sum_{\langle i,j \rangle} s_i s_j,
$$

where

- $s_i = \pm 1$ is the spin at lattice site $i$,
- $J$ is the nearest-neighbor coupling,
- the sum runs over nearest-neighbor pairs.

The simulations use periodic boundary conditions, so that the square lattice has no physical edges.

The calculations use units in which

$$
k_B = 1.
$$

For a ferromagnetic system, $J>0$, neighboring spins energetically favor alignment.

The exact critical temperature of the infinite square-lattice Ising model is

$$
T_c =
\frac{2J}
{\ln(1+\sqrt{2})}.
$$

For $J=1$,

$$
T_c \approx 2.269.
$$

### Metropolis dynamics

The system is evolved using the Metropolis Monte Carlo algorithm.

At every attempted update, one lattice site is selected randomly and its spin is proposed to flip,

$$
s_i \rightarrow -s_i.
$$

Since only the interaction of that spin with its four nearest neighbors changes, the corresponding energy variation is

$$
\Delta E
=
2J s_i
\sum_{j \in \mathrm{nn}(i)} s_j.
$$

The proposed spin flip is accepted with probability

$$
P_{\mathrm{accept}} =
\begin{cases}
1, & \Delta E \leq 0, \\
e^{-\Delta E/T}, & \Delta E > 0.
\end{cases}
$$

Energy-lowering moves are therefore always accepted, while energy-increasing moves can still occur because of thermal fluctuations.

At low temperature, unfavorable moves are strongly suppressed and the system tends toward ordered configurations. At high temperature, thermal fluctuations make spin flips much more frequent and destroy long-range magnetic order.

One Monte Carlo cycle in this project corresponds to

$$
L^2
$$

attempted spin flips, where $L$ is the linear lattice size. Therefore one cycle corresponds, on average, to one attempted update per lattice site.

## Results

The project can be used either to study the time evolution of a single system at fixed temperature or to investigate thermodynamic behavior over a range of temperatures.

### Single-temperature simulation

A single simulation records the magnetization during the Monte Carlo evolution.

![Magnetization as a function of Monte Carlo cycles](figures/single_magnetization.png)

The final spin configuration can also be visualized directly.

![Final Ising lattice](figures/single_final_lattice.png)

The time evolution makes it possible to observe the equilibration process and the fluctuations around the equilibrium state.

### Temperature sweep

Multiple independent simulations can be performed over a range of temperatures.

The main observable is the absolute magnetization per site,

$$
|m|
=
\left|
\frac{1}{N}
\sum_i s_i
\right|.
$$

The Monte Carlo results are compared with the exact spontaneous magnetization of the infinite two-dimensional Ising model.

![Magnetization versus temperature](figures/magnetization_vs_temperature.png)

The transition from the ordered phase to the disordered phase occurs around the exact critical temperature $T_c$.

For a finite lattice, the measured value of $\langle |m| \rangle$ does not become exactly zero above $T_c$. This is expected because the simulation contains a finite number of spins.

Even in the disordered phase, instantaneous fluctuations generally produce a small non-zero total magnetization. Since the observable used here is the absolute magnetization, positive and negative fluctuations do not cancel each other.

Therefore,

$$
\langle |m| \rangle > 0
$$

for a finite system even above the critical temperature.

The exact theoretical curve instead describes the thermodynamic limit

$$
L \rightarrow \infty,
$$

where the spontaneous magnetization is exactly zero for $T \geq T_c$.

The energy per site is also calculated from the final configurations and compared with the exact thermodynamic-limit result.

![Energy versus temperature](figures/energy_vs_temperature.png)

In the current implementation, the energy is evaluated only from the final configuration of each independent simulation. It is therefore displayed without a statistical error estimate.

## Project structure

```text
.
├── configs/
│   ├── single_run_parameters.toml
│   └── multi_run_parameters.toml
├── figures/
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