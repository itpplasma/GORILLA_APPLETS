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
Current branch `self_consistent_electric_field` focuses on:
- Electron diffusion coefficient computation (`flux_deviation_mod.f90`)
- Self-consistent electric field calculations (`utils_self_consistent_ef_mod.f90`)
- Enhanced particle collision modeling
- Flux surface analysis

### Data Visualization
- **MATLAB/**: Comprehensive plotting scripts for each computational mode
- Output format: ASCII data files, NetCDF for field data
- Example: `mono_energetic_transport.m` for transport coefficient visualization

### Testing Approach
- Integration tests via EXAMPLES/ directory
- Physics validation: reversibility, Boltzmann distribution, Poincaré invariant
- No formal unit test framework; use `make test` for available tests