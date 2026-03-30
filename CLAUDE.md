# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview
GORILLA_APPLETS is a scientific computing application for plasma physics that computes guiding-center orbits for charged particles in toroidal fusion devices with 3D field geometry. It extends the GORILLA (Guiding-center ORbit Integration with Local Linearization Approach) code with specialized applications.

## Build Commands
```bash
# Standard build (using CMake with Ninja)
make

# Alternative build configurations
make CONFIG=Debug      # Debug build with extra checks
make CONFIG=Release    # Optimized release build (default)

# Clean build
make clean

# Run tests (if available)
make test

# Reconfigure CMake
make reconfigure
```

The build produces `gorilla_applets_main.x` executable which supports 12 different computational modes selected via `i_option` parameter.

## Architecture Overview

### Core Structure
- **SRC/**: Main application source code (Fortran 90)
  - `gorilla_applets_main.f90`: Entry point, dispatches to different computational modes
  - `gorilla_applets_settings_mod.f90`: Global configuration management
  - `gorilla_applets_sub_mod.f90`: Core computational routines
  - Various `*_mod.f90` files: Specialized modules for each computational mode
- **SRC_CORE**: Symlink to GORILLA core functionality (../GORILLA/SRC)
- **BUILD/**: CMake build directory (generated)
- **INPUT/**: Configuration files for different applications
- **EXAMPLES/**: Ready-to-run test cases with configurations

### Computational Modes (i_option)
1. Pre-computation of fluxtube volume
2. Mono-energetic transport coefficient calculation
3. Collisionality scan for transport coefficients
4. Numerical diffusion coefficient
5. Alpha particle lifetime calculation
6. Direct VMEC integrator
7. Poincaré invariance computation
8. Reversibility test
9. Boltzmann test
10. Field line tracing
11. Divertor heat loads calculation
12. Self-consistent electric field computation

### Key Dependencies
- **Fortran compiler** with OpenMP support
- **BLAS/LAPACK**: Linear algebra operations
- **NetCDF**: Scientific data I/O for equilibrium files
- **CMake >= 3.12**: Build system
- **Ninja**: Build generator (preferred)

### Running the Code
```bash
# Basic execution
./gorilla_applets_main.x

# Control parallelization
export OMP_NUM_THREADS=<number_of_cores>  # Avoid hyperthreading
./gorilla_applets_main.x
```

### Required Input Files
Place in working directory:
- `gorilla_applets.inp`: Main configuration (selects i_option)
- `tetra_grid.inp`: Tetrahedral grid settings
- `gorilla.inp`: Core GORILLA settings
- Mode-specific files from INPUT/ directory based on i_option
- MHD equilibrium files (NetCDF or g-file format)

### Active Development Areas

**Branch `transport-km-benchmark`**: Kramers-Moyal displacement-based D11 benchmark
- Measures D11(s) via `<(ds)^2>/(2dt)` per flux surface (Kramers-Moyal / Einstein relation)
- Existing infrastructure: `calc_electron_diffusion_coefficients` at `utils_self_consistent_ef_mod.f90:1254-1422`
- KM fitting: `fit_transport_coefficients` at `transport_benchmark_utils_mod.f90:22-64`
- Benchmark target: GORILLA vs NEO-2-QL vs NEOART on AUG 39461 neon discharge

**Archived branch `global-transport-fitting`**: LM inverse problem approach (abandoned)
- Steady-state PDE fitting failed: MC data doesn't reach steady state in feasible tracing time
- Conservation violated by 100-30000x (net flux out << total source weight)
- The displacement-based approach avoids this entirely

### Transport Coefficient Estimation: Kramers-Moyal Method

For a diffusion process ds = sqrt(2D*dt)*xi + A*dt:
- **D11(s) = <(ds)^2> / (2dt)** at each flux surface (2nd Kramers-Moyal coefficient)
- **A(s) = <ds> / dt** at each flux surface (1st coefficient = drift/convection)
- No PDE solving, no iterative fitting, O(N) computation
- Works with any tracing time in the diffusive regime (no steady-state needed)

Existing code accumulates `s%delta_s` and `s%delta_s_squared` during particle pushing
(when `boole_diffusion_coefficient=.true.`), then `fit_transport_coefficients` extracts
A and D11 from the linear/quadratic growth of these moments vs time.

### Benchmark: GORILLA vs NEO-2-QL vs NEOART

Existing NEO-2-QL vs NEOART benchmark at AUG shots 39084/39461 (neon discharges):
- Data: `$DATA/AUG/NEO-2/39461/neoart_benchmark_neon_discharge/`
- NEOART data: `$DATA/AUG/NEO-ART/39461/neoart_benchmark_neon_discharge/`
- Comparison tools: `$DATA/AUG/NEO-2/neo2_neoart_benchmark_tools/`
- Conversion: D11_AX = D11_s / <nabla_s>^2 (see benchmark tools README)

### Data Visualization
- **MATLAB/**: Legacy plotting scripts for each computational mode
- Output format: ASCII data files, NetCDF for field data

### Testing
- `make test` runs unit tests
- `test_global_transport_fit`: LM fit tests (manufactured profiles, adjoint, boundary conditions)
- `test_kramers_moyal`: KM displacement tests (1D random walk recovery, fit_transport_coefficients)