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
| 13 | Anomalous transport |
| 14 | Helical core particle tracing |

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


### CINECA Leonardo

Built and tested on the Leonardo CPU login nodes and in a budget-free serial
batch job (GCC 12.2.0, CMake 3.27.9, OpenBLAS 0.3.26, NetCDF-C 4.9.2,
NetCDF-Fortran 4.6.1, Python 3.11.7) at GORILLA_APPLETS `cb6f347` with GORILLA
`1a6cdc1`, 2026-10-19. Pitagora access was not available, so nothing is
verified there. GitHub is reachable from the login nodes. The shims below are
needed because the Spack NetCDF modules ship no `netCDFConfig.cmake` and the
CMake files link `-lblas` and `-llapack`, which OpenBLAS does not provide as
separate libraries. The Python tests need `f90nml`, installed into a virtual
environment.

```bash
module purge
module load profile/base gcc/12.2.0 cmake/3.27.9 ninja/1.11.1 \
    openblas/0.3.26--gcc--12.2.0 python/3.11.7 \
    netcdf-c/4.9.2--gcc--12.2.0-spack0.22 \
    netcdf-fortran/4.6.1--gcc--12.2.0-spack0.22
export B=$CINECA_SCRATCH/builds TMPDIR=$CINECA_SCRATCH/builds/tmp
mkdir -p $TMPDIR $B/shim/netCDF $B/shim/lib
export CMAKE_BUILD_PARALLEL_LEVEL=4 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
echo 'set(netCDF_FOUND TRUE)' > $B/shim/netCDF/netCDFConfig.cmake
for l in lapack blas; do ln -sf $OPENBLAS_LIB/libopenblas.so $B/shim/lib/lib$l.so; done
LINK="-L$B/shim/lib -Wl,-rpath,$OPENBLAS_LIB -L$NETCDF_FORTRAN_LIB -L$NETCDF_C_LIB -Wl,-rpath,$NETCDF_FORTRAN_LIB -Wl,-rpath,$NETCDF_C_LIB"
INC="-I$NETCDF_FORTRAN_INC -I$NETCDF_C_INC"
cd $B
git clone https://github.com/itpplasma/GORILLA.git
git clone https://github.com/itpplasma/GORILLA_APPLETS.git && cd GORILLA_APPLETS
python3 -m venv $B/venv && $B/venv/bin/pip install f90nml
nice cmake -S . -B build-lrd -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DPython3_EXECUTABLE=$B/venv/bin/python \
    -DnetCDF_DIR=$B/shim/netCDF "-DCMAKE_Fortran_FLAGS=$INC" \
    "-DCMAKE_EXE_LINKER_FLAGS=$LINK" "-DCMAKE_SHARED_LINKER_FLAGS=$LINK"
nice cmake --build build-lrd -j 4
(cd build-lrd && ctest --output-on-failure)
```

The sibling `GORILLA` checkout is required. Release build: 36 s with `-j 4`
under `nice`. All 9 registered tests passed (32 s in total, four threads;
the three Fortran unit tests take under a second, the Python application
tests up to 16 s each). Run `ctest` in a batch job on `--partition=lrd_all_serial`
(budget-free, no account, `--cpus-per-task=4`) rather than on the login node.

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

## Merging

- Pull requests, including every layer of a stacked PR, run the fast CI tier.
  Label a PR `full-ci` to run the full tier on it before queueing.
- Review a stack layer by layer. When it is ready, collapse it into one PR into
  `main` (retarget the top layer to `main`); the lower PRs then show as merged.
- Click "Merge when ready" (or `gh pr merge --merge --auto`). The PR enters the
  merge queue, the full CI tier runs on the exact merge result, and on success
  it lands on `main` as a merge commit. `main` itself is not re-tested on push.
- A branch need not be up to date with `main`. Merge `main` into a PR only
  when it conflicts.

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
