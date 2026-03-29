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
Current branch `global-transport-fitting` focuses on:
- Global transport coefficient fitting via Levenberg-Marquardt inverse model
- Local transport reference measurement via delta-f Monte Carlo
- Internal validation: global vs local A(s), B(s) comparison with 2-sigma bands
- External benchmark preparation: GORILLA vs NEO-2-QL/NEO-2-PAR/NEOART

### Global Transport Fitting (i_option=15)

The global transport fit solves an inverse problem to recover radial transport coefficients A(s) and B(s) from Monte Carlo particle-tracing experiments on the full domain.

#### Physics: Neoclassical Transport Coefficients

The neoclassical particle flux in the linearized drift-kinetic framework (Schatzlmayr thesis, Eq. 285-288) is:

```
Gamma = -D11 * A1 - D12 * A2
```

where the thermodynamic forces are:

```
A1 = grad(n)/n - 3/2 * grad(T)/T - q*E/T    (density + E-field drive)
A2 = grad(T)/T                                 (temperature drive)
```

and D11, D12 are defined through the linearized drift-kinetic operator:

```
n * D11 * e_j = integral d3v v * L_hat^{-1}(v . e_j * f0)
n * D12 * e_j = integral d3v v * L_hat^{-1}(v . e_j * m*v^2/(2T) * f0)
```

In the mono-energetic Lorentz operator limit (no T gradient, no E field):
- A1 = d(ln n)/ds,  A2 = 0
- Gamma/n = -D11 * d(ln n)/ds
- Comparing with our model Gamma/n = A - B * d(ln n)/ds: **B = D11**, **A = 0**

#### Code Structure

Key modules for the global transport fit:

| Module | Purpose |
|--------|---------|
| `global_transport_fit_gorilla_mod.f90` | Orchestration: MC experiments, local reference, output |
| `global_transport_fit_core_mod.f90` | LM inverse model, forward problem, adjoint gradient |
| `global_transport_fit_types_mod.f90` | Data types for experiments, control, results |
| `global_transport_fit_settings_mod.f90` | Namelist-driven configuration |
| `global_transport_fit_io_mod.f90` | File I/O for experiments, profiles, convergence history |
| `global_transport_fit_math_mod.f90` | Geometry, basis functions, linear algebra |
| `transport_statistics_mod.f90` | Signal quality, local fit statistics, sanitization |
| `usual_transport_benchmark_mod.f90` | Local MC transport measurement at single surface |

#### Delta-f Particle Weighting Method

The density gradient in local transport measurements is encoded as a **delta-f particle weight**, not a physical change in the background plasma:

```
weight = exp(gradient * (s - s_ref))
```

Applied in `calc_particle_weights_and_jperp()` at `utils_self_consistent_ef_mod.f90`.

The collision operator uses a separate, fixed background density (`background_density_cm3`). The collision frequency `nu ~ L^{alpha,beta} * n_beta / v_th^3` depends on the background density, NOT the particle weight. This is the standard approach used by NEO-2, FORTEC-3D, and SIMPLE.

For each local measurement at surface s*:
- `density_profile_reference_s = s*` (reference IS the measurement point)
- `background_density_cm3 = in%density` (fixed, same as global experiment)
- The gradient only enters through the particle weight
- `normalized_flux = weighted_flux_density / in%density` gives Gamma/n(s*)

#### Forward Problem

The discrete 1D transport equation on flux surfaces:

```
div(Phi) = V * S    (mass conservation per shell)
Phi(j) = A(s_j) * 0.5*(n_{j-1} + n_j) - B(s_j) * (n_j - n_{j-1}) / delta_s
```

Boundary conditions:
- Inner (magnetic axis): Phi(1) = 0 (zero area, toroidal topology)
- Outer: prescribed density (absorbing boundary)

B(s) is parameterized as exp(linear basis) to ensure positivity.

#### Collision Operators

From Schatzlmayr thesis Chapter 4, the Fokker-Planck collision operator is implemented as a diffusion process in velocity space with:

| Operator | Code value | Description |
|----------|-----------|-------------|
| Lorentz | 4 | Pitch-angle scattering only (mono-energetic), `lag=0` internally |
| Full Fokker-Planck | other | Energy + pitch diffusion + drag |

The Lorentz operator (collision_operator=4) is mono-energetic: particles keep their energy, only the pitch angle lambda changes. This is the simplest case for benchmarking D11.

Collision coefficients (Eq. 251, 255 in thesis):
```
D_{lambda,lambda} = nu_0^{alpha,beta} * (1-lambda^2)/v_tilde^2 * A
V_lambda = -2*lambda/(1-lambda^2) * D_{lambda,lambda}
```
where A, B, C are functions of u = v_alpha / v_{beta,th} using erf(u) and erf'(u).

#### Signal Quality Gating

Both source experiments must pass quality gates before fitting:
- `density_source_relstd >= 5%` (density and source profiles must differ)
- `n_supported_flux_boundaries >= 3` (enough radial flux signal above 2-sigma noise)
- Retry with increased tracing time if not met (up to max_source_trials)

### Reference Documents

| Document | Location | Content |
|----------|----------|---------|
| Schatzlmayr PhD thesis | `~/Nextcloud/tug/stud/PhD_Jonatan_Schatzlmayr/dissertation_jonatan_schatzlmayr.pdf` | GORILLA theory, symplecticity, collision operator, self-consistent E-field |
| Eder PhD thesis | `~/Nextcloud/plasma/DOCUMENTS/Eder_Dissertation/eder_diss_final.pdf` | Original GORILLA geometric integration, orbit benchmarks |
| GORILLA PoP 2020 | `~/Nextcloud/plasma/DOCUMENTS/GORILLA/2020_PoP_ver03/gorilla_paper.pdf` | Published GORILLA paper with transport benchmarks |
| NEO-2 notes | `~/Nextcloud/plasma/DOCUMENTS/INFO/codes/NEO-2/NEO-2.md` | NEO-2 input/output documentation |
| Transport benchmark data | `~/Nextcloud/plasma/DOCUMENTS/GORILLA/EDER_2020/Plots/RAW_Data/transport_coefficient/` | Existing D11 collisionality/E-field scans |
| 2025 benchmark suite | `~/Nextcloud/plasma/DOCUMENTS/eu_us_transport_task_force_workshop/bootstrap_ttf_2025/graphs/benchmark/` | Recent axisymmetric/QH transport benchmarks |
| DKES comparison data | `~/Nextcloud/plasma/DOCUMENTS/NTV_Notes_Andreas/NTV_Paper_Mar2014/DRAFT_SECIV/DATA/DKES_NEW/` | Historical D11/D12 benchmark data |

### Data Visualization
- **MATLAB/**: Legacy plotting scripts for each computational mode
- **scripts/**: Python validation and plotting (untracked, local-only)
  - `run_global_transport_validation.sh` - end-to-end validation workflow
  - `plot_global_transport_validation.py` - generates A(s), B(s), signal quality, convergence plots
- Output format: ASCII data files, NetCDF for field data

### Testing
- `make test` runs the unit test suite (`TESTS/test_global_transport_fit.f90`)
- 13 tests covering: manufactured profile recovery, adjoint gradient, axis BC, geometry, signal quality, sanitization, local flux pair recovery, boundary selection, convergence history
- All tests must pass 100% before commits