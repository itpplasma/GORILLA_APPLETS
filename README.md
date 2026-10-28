# GORILLA_APPLETS

Plasma-physics applications built on the guiding-center orbit integrator
[GORILLA](https://github.com/itpplasma/GORILLA) (**G**uiding-center **OR**bit
**I**ntegration with **L**ocal **L**inearization **A**pproach).

GORILLA traces guiding-center orbits in three-dimensional toroidal fields with
a quasi-geometric integrator that conserves total energy, magnetic moment and
phase-space volume (Refs. [1, 2]). GORILLA_APPLETS provides a single executable,
`gorilla_applets_main.x`, that runs one of several applications on top of it.

## Applications

The application is selected with `i_option` in the namelist
`gorilla_applets_nml` of `gorilla_applets.inp`:

1. Pre-computation of flux-tube volume
2. Mono-energetic transport coefficient for one collisionality
3. Collisionality scan of mono-energetic transport coefficients
4. Numerical diffusion coefficient (no collisions)
5. Alpha-particle lifetime
6. Direct VMEC integrator
7. Poincaré invariance test
8. Reversibility test
9. Boltzmann test
10. Field-line tracing
11. Divertor heat loads
12. Self-consistent electric field
13. Anomalous transport
14. Helical-core particle tracing
15. RMP response currents

## Building

The build expects a checkout of GORILLA next to this repository; `SRC_CORE`
is a symbolic link to `../GORILLA/SRC`.

```bash
git clone https://github.com/itpplasma/GORILLA.git
git clone https://github.com/itpplasma/GORILLA_APPLETS.git
cd GORILLA_APPLETS
make                  # Release build with CMake and Ninja
```

Requirements: a Fortran compiler with OpenMP, BLAS/LAPACK, NetCDF with
Fortran bindings, CMake 3.12 or newer and Ninja. The tests also need Python 3
with `f90nml`. The executable is `BUILD/gorilla_applets_main.x`.
See [docs/building.md](docs/building.md) for build configurations, tests and
notes for the CINECA Leonardo cluster.

## Running

Run the executable in a directory that holds the input files. Every folder in
[EXAMPLES/](EXAMPLES/) links to `BUILD/gorilla_applets_main.x`:

```bash
cd EXAMPLES/ALPHA_LIFETIME
export OMP_NUM_THREADS=<physical cores>
./gorilla_applets_main.x
```

The working directory needs:

- `gorilla.inp` and `tetra_grid.inp`: GORILLA settings (templates in
  `GORILLA/INPUT/`),
- `gorilla_applets.inp`: selects `i_option`,
- the input file of the selected application, for example
  `alpha_lifetime.inp` or `mono_energetic_transp_coef.inp`,
- the MHD equilibrium (VMEC NetCDF or g-file) and, for applications that use
  random numbers, `seed.inp`.

Templates are in [INPUT/](INPUT/). The folders in [EXAMPLES/](EXAMPLES/) contain
complete setups for most applications, including input files, links to
equilibria and MATLAB or Python plotting scripts.

## Testing

```bash
make test
```

runs the CTest suite defined in [TESTS/](TESTS/) and `CMakeLists.txt`: Fortran
unit tests and short application runs (small grids, few particles) driven by
Python scripts. The [CI workflow](.github/workflows/CI.yml) builds against a
fresh checkout of GORILLA and runs the suite on pull requests.

## Documentation

- [docs/building.md](docs/building.md): build configurations, tests, CINECA
  Leonardo.
- [docs/contributing.md](docs/contributing.md): CI tiers and merging.
- [docs/plans/](docs/plans/): design notes for the Ornstein-Uhlenbeck
  collision operator.
- For the GORILLA core (method, input files, `orbit_timestep_gorilla`,
  `gorilla_plot`), see the [GORILLA repository](https://github.com/itpplasma/GORILLA).

## License

MIT License, see [LICENSE](LICENSE). GORILLA is also MIT-licensed and bundles a
few third-party routines under their own licences; see its README.

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
