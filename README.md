# GORILLA_APPLETS

Applications built on top of [GORILLA](https://github.com/itpplasma/GORILLA) (**G**uiding-center **OR**bit **I**ntegration with **L**ocal **L**inearization **A**pproach).

GORILLA itself traces guiding-center orbits in three-dimensional toroidal fields via a quasi-geometric integrator that conserves total energy, magnetic moment, and phase-space volume; see Refs. [1,2] and the [GORILLA repository](https://github.com/itpplasma/GORILLA) for the underlying method and core interfaces (`orbit_timestep_gorilla`, `gorilla_plot`).

GORILLA_APPLETS wraps the GORILLA core with a set of higher-level plasma-physics applications: a single executable `gorilla_applets_main.x` dispatches to one of several computational modes (alpha-particle losses, mono-energetic transport, Poincaré-invariance and reversibility diagnostics, divertor heat loads, self-consistent electric field, and more) selected via the `i_option` key in `gorilla_applets.inp`.

## Applications

Selected at runtime via `gorilla_applets_nml/i_option` in `gorilla_applets.inp`:

| `i_option` | Application |
|------------|-------------|
| 1  | Pre-computation of fluxtube volume |
| 2  | Mono-energetic transport coefficient |
| 3  | Collisionality scan for transport coefficients |
| 4  | Numerical diffusion coefficient |
| 5  | Alpha-particle lifetime |
| 6  | Direct VMEC integrator |
| 7  | Poincaré-invariance test |
| 8  | Reversibility test |
| 9  | Boltzmann test |
| 10 | Field-line tracing |
| 11 | Divertor heat loads |
| 12 | Self-consistent electric field |

## Building

GORILLA_APPLETS expects a sibling checkout of GORILLA next to this repository (the build symlinks `SRC_CORE` to `../GORILLA/SRC`):

```bash
git clone git@github.com:itpplasma/GORILLA.git
git clone git@github.com:itpplasma/GORILLA_APPLETS.git
cd GORILLA_APPLETS
make                  # Release build via CMake + Ninja
make CONFIG=Debug     # Debug build with extra checks
make clean
make reconfigure
```

The resulting executable is `BUILD/gorilla_applets_main.x` (symlinked to the project root). Requires a Fortran compiler with OpenMP, BLAS/LAPACK, NetCDF (Fortran bindings), CMake >= 3.12, and Ninja.

## Running

```bash
export OMP_NUM_THREADS=<physical_cores>   # avoid hyperthreading
./gorilla_applets_main.x
```

The working directory needs the following inputs (blueprints under [INPUT/](INPUT/) and in the sibling GORILLA repo):

- `gorilla.inp`, `tetra_grid.inp` — core GORILLA settings (from GORILLA/INPUT/)
- `gorilla_applets.inp` — selects `i_option` and global applet settings
- One mode-specific input file (e.g. `alpha_lifetime.inp`, `mono_energetic_transp_coef.inp`, `boltzmann.inp`, `reversibility_test.inp`, `divertor_heat_loads.inp`, `field_line_tracing.inp`, `self_consistent_ef.inp`)
- MHD equilibrium files (VMEC NetCDF or g-file), and a `seed.inp` where the application uses random numbers

The runnable scenarios in [EXAMPLES/](EXAMPLES/) show the full setup per application; each folder contains its input files, MHD equilibrium symlinks, and a MATLAB or Python plotting script for the output.

## Testing

Continuous integration runs out of [TESTS/](TESTS/) via CTest:

```bash
make test
```

Each smoke test is a short scenario (small grid, few particles) wired through CMake; the alpha-particle lifetime test, for example, drives the binary twice (flux-tube precomputation, then a 10-particle trace) and checks the output file. The GitHub Actions workflow in [.github/workflows/CI.yml](.github/workflows/CI.yml) builds against a fresh checkout of GORILLA and runs the suite on every PR.

## License

GORILLA_APPLETS is released under the MIT License — see [LICENSE](LICENSE). The build depends on the GORILLA core ([itpplasma/GORILLA](https://github.com/itpplasma/GORILLA)), which is also MIT-licensed and itself bundles a small number of third-party routines under their own licenses; see GORILLA's README for details.

## References

When using this code in scientific publications, please cite [1] and [2]:

[1] M. Eder, C.G. Albert, L.M.P. Bauer, S.V. Kasilov and W. Kernbichler,
*Quasi-geometric integration of guiding-center orbits in piecewise linear toroidal fields*,
Physics of Plasmas **27**, 122508 (2020).
<https://doi.org/10.1063/5.0022117> — Preprint: <https://arxiv.org/abs/2007.08151>

[2] M. Eder, C.G. Albert, L.M.P. Bauer, G.S. Graßler, S.V. Kasilov, W. Kernbichler, M. Meisterhofer and M. Scheidt,
*GORILLA: Guiding-center ORbit Integration with Local Linearization Approach*,
Journal of Open Source Software.
Preprint: <https://github.com/openjournals/joss-papers/blob/joss.03116/joss.03116/10.21105.joss.03116.pdf>
