# OpenQuantumSpinDynamics

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
├── Manifest.toml                 # Dependencies for the project
├── Project.toml                 # Project metadata and dependencies
├── src                          # Source code
│   ├── Arnoldi.jl              # Implementation of the Arnoldi algorithm
│   ├── Coupling.jl             # Tools for handling system coupling
│   ├── Disorder.jl             # Modules for introducing disorder in systems
│   ├── KrylovArnoldiSolver.jl  # Efficient Krylov space solver
│   ├── LindbladSolver.jl       # Solver for Lindblad master equations
│   ├── OpenQuantumSpinDynamics.jl # Main library module
│   ├── PauliOps.jl             # Pauli operators utilities
│   ├── QuantumState.jl         # Representation and operations on quantum states
│   ├── Solvers.jl              # General solver implementations
│   ├── SpinModels.jl           # Models for quantum spin systems
│   ├── StochasticWavefunctionSolver.jl # Solver for stochastic wavefunctions
│   ├── setup.jl                # Setup configurations
│   └── utils.jl                # Helper functions and utilities
├── test                        # Test suite
│   ├── configuration_test.json # Sample configuration for tests
│   └── runtests.jl             # Test runner
├── main.jl                     # Main entry point for running simulations
└── plot.jl                     # Script for generating plots
```

## Quick Start

### Prerequisites

Ensure that you have Julia installed on your system. The recommended version is Julia 1.8 or higher.

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
julia --project=../OpenQuantumSpinDynamics -t 4 main.jl
```

### Generating Plots

Generate visualization of the simulation results:
```bash
julia --project=../OpenQuantumSpinDynamics plot.jl
```

## Contribution Guidelines

Contributions are welcome! To contribute:
1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Submit a pull request with a clear description of your changes.

Please ensure all changes pass the tests before submitting:
```julia
julia --project=../OpenQuantumSpinDynamics test/runtests.jl
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

Special thanks to the contributors and the Julia community for their support and inspiration.


