# OpenQuantumSpinDynamics

[![CI Pipeline](https://github.com/javahedi/OpenQuantumSpinDynamics/actions/workflows/ci.yml/badge.svg)](https://github.com/javahedi/OpenQuantumSpinDynamics/actions/workflows/ci.yml)

![Code Style](https://img.shields.io/badge/code%20style-JuliaFormatter-blue)

OpenQuantumSpinDynamics is a Julia-based library designed for simulating the dynamics of open quantum spin systems. The library provides robust tools for solving Lindblad master equations, stochastic wavefunctions, and analyzing spin models with a variety of solvers. 

This project leverages the power of Julia's ecosystem for high-performance numerical computing and is optimized for parallel execution.

## Features

- **Comprehensive Spin Dynamics Tools**: 
  - Simulation of quantum spin systems.
  - Lindblad master equation solver.
  - Stochastic wavefunction solver.

- **Powerful Solvers**: 
  - Krylov Arnoldi solver for efficient computation.
  - Utilities for handling coupling, disorder, and Pauli operators.

- **Parallel Computing Support**: 
  - Distributed and parallel processing using Julia's `Distributed` and `pmap`.
  - GPU acceleration for Monte Carlo sampling in `StochasticWavefunctionSolver.jl`.

- **Customization and Flexibility**: 
  - Modular design to adapt to various research needs.
  - Easily configurable setup and test structure.

## Repository Structure

```
OpenQuantumSpinDynamics
    ├── Manifest.toml
    ├── Project.toml
    ├── README.md
    ├── examples
    │   ├── configuration.json
    │   ├── main.jl
    │   └── plot.jl
    ├── src
    │   ├── Arnoldi.jl
    │   ├── Coupling.jl
    │   ├── Disorder.jl
    │   ├── KrylovArnoldiSolver.jl
    │   ├── LindbladSolver.jl
    │   ├── OpenQuantumSpinDynamics.jl
    │   ├── PauliOps.jl
    │   ├── QuantumState.jl
    │   ├── Solvers.jl
    │   ├── SpinModels.jl
    │   ├── StochasticWavefunctionSolver.jl
    │   ├── setup.jl
    │   └── utils.jl
    └── test
        ├── configuration_test.json
        └── runtests.jl
```

## Quick Start

### Prerequisites

Ensure that you have Julia installed on your system. The recommended version is Julia 1.11 or higher.

### Installation

Clone the repository:
```bash
git clone https://github.com/javahedi/OpenQuantumSpinDynamics.git
cd OpenQuantumSpinDynamics
```

Activate the project environment in Julia:
```julia
julia --project=.
using Pkg
Pkg.instantiate()
```

This will install all required dependencies listed in the `Project.toml` file.

### Running Simulations

Run the main simulation script with 4 threads:
```bash
julia --project=. -t 4 examples/main.jl   
```

### Generating Plots

Generate visualization of the simulation results:
```bash
julia --project=. examples/plot.jl   
```

## Contribution Guidelines

Contributions are welcome! To contribute:
1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Submit a pull request with a clear description of your changes.

Please ensure all changes pass the tests before submitting:
```julia
julia --project=. test/runtests.jl
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

Special thanks to the contributors and the Julia community for their support and inspiration.


