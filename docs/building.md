# Building GORILLA_APPLETS

## Layout

GORILLA_APPLETS compiles the GORILLA core sources together with its own. The
repository contains a symbolic link `SRC_CORE -> ../GORILLA/SRC`, so a GORILLA
checkout must sit next to this repository:

```bash
git clone https://github.com/itpplasma/GORILLA.git
git clone https://github.com/itpplasma/GORILLA_APPLETS.git
```

## Requirements

- Fortran compiler with OpenMP (GNU Fortran is used in CI)
- BLAS and LAPACK
- NetCDF-C and NetCDF-Fortran
- CMake 3.12 or newer and Ninja
- Python 3 with `f90nml` for the application tests

On macOS, `CMakeLists.txt` locates NetCDF through `nf-config` and `nc-config`,
so both Homebrew and MacPorts installations work.

## Make targets

The top-level `Makefile` wraps CMake with the Ninja generator in `BUILD/`:

```bash
make                     # Release build (default CONFIG=Release)
make CONFIG=Debug        # Debug build with extra checks
make reconfigure CONFIG=Debug   # change the build type of an existing BUILD/
make test                # ctest --output-on-failure in BUILD/
make clean               # remove BUILD/
```

The executable is `BUILD/gorilla_applets_main.x`. Every example folder links
to it directly.

## Batched transport with CUDA

Applications selected by `i_option=2`, `3` or `4` in
`gorilla_applets_nml` can advance their active markers in batches. Set
`transport_backend` in the `TRANSPCOEFNML` namelist of
`mono_energetic_transp_coef.inp` to `cpu` (the default scalar tracer),
`cpu-batch` (the OpenMP batch reference) or `cuda`.

Build against the GORILLA revision recorded in `GORILLA_REF`. CUDA is optional
and needs CMake 3.18 or newer, a CUDA toolkit and a supported NVIDIA GPU:

```bash
git -C ../GORILLA checkout "$(cat GORILLA_REF)"
cmake -S . -B BUILD -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DGORILLA_ENABLE_CUDA=ON -DGORILLA_RUNTIME_CHECKS=OFF \
  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_CUDA_COMPILER=nvcc \
  -DCMAKE_CUDA_ARCHITECTURES=75
cmake --build BUILD
OMP_NUM_THREADS=4 ctest --test-dir BUILD --output-on-failure \
  -R 'mono_energetic_transport_(slow|cuda)'
```

Replace architecture `75` (Tesla T4) with the target GPU's architecture,
for example `120` for Blackwell with a toolkit that supports it. A CPU-only
build needs no CUDA installation.

The CUDA path currently requires polynomial order 2, `i_precomp=0`, constant
time parametrization (`i_time_tracing_option=1`),
`handover_processing_kind=1`, `boole_adaptive_time_steps=.false.` and
`boole_strong_electric_field=.false.`. Disable optional orbit quantities with
`boole_time_Hamiltonian=.false.`, `boole_gyrophase=.false.`,
`boole_vpar_int=.false.` and `boole_vpar2_int=.false.`. Batch transport requires
`i_integrator_type=1`, `boole_random_precalc=.true.` and no helical perturbation.
Unsupported settings fail explicitly.

Lorentz scattering, measurements and transport estimation stay on the host.
A cell rejected by the CUDA checks is retraced from its unchanged entry state
by the CPU pusher. Numerical failures abort the estimate; physical losses
retain all successful earlier samples. Timings report setup, transfers,
device execution, CPU recovery and complete batch-consumer wall time.

## Tests

The CTest suite consists of Fortran unit tests (collision operators, RMP
loading and profile moments, transport normalisation and loss moments) and
short application runs driven by the Python scripts in `TESTS/`, for example
the alpha-particle lifetime test, which runs the flux-tube volume
precomputation and then a ten-particle trace and checks the output.
`TESTS/OU_COLLISION_INVARIANT.md` and `TESTS/RMP_DENSITY_RESPONSE.md` describe
two of the physics checks.

## Transport benchmarks

`TESTS/benchmark_transport.py` measures the actual `cpu`, `cpu-batch` and
`cuda` transport backends. Build with CUDA enabled to include `cuda`, and run
the native public-equilibrium test first to prepare its input fixture:

```bash
ctest --test-dir BUILD -R '^mono_energetic_transport_slow$' --output-on-failure
python3 TESTS/benchmark_transport.py \
    --binary BUILD/gorilla_applets_main.x \
    --fixture BUILD/TESTS/mono_energetic_transport_work \
    --output transport-benchmark \
    --particles 128 1024 8192 --threads 1 2 4 8 \
    --backends cpu cpu-batch cuda --steps 32 --time-step 1e-9 \
    --species electron --repeats 3 --plot
```

The output directory must be new. Each process gets a separate copy of the
fixture with the same seed, second-order polynomial pusher, energy, duration
and diagnostics. Collisionless option 4 uses the requested number of
measurements and spacing in seconds. For seeded Lorentz collisions, use
`--collisions lorentz --flight-time-multiplier 1` without `--steps` or
`--time-step`: option 2 resolves its cadence and duration from collision and
bounce times. Electron and deuterium runs are separate benchmark invocations.
To distinguish changes in CPU mathematics from GPU acceleration, add
`--baseline-binary /path/to/older/gorilla_applets_main.x`. Its scalar runs are
labelled `baseline-cpu`; `cpu` always uses the current `--binary`. The older
input omits the new backend selector. Record its revision with
`--source-revision baseline=REV` and its compiler with `--baseline-compiler`.

JSON and CSV retain every process, its log, terminal counters and available
setup, upload, kernel, download and recovery timings. Whole-process wall time
includes mesh setup and diagnostics. The first process for each configuration
is recorded separately; later repeats are fresh processes that may benefit
from filesystem caches, rather than a resident GPU warmup. Medians and ranges
use the later repeats. Output comparisons use the transport test's tolerances;
failed processes, non-finite output or mismatched observables invalidate the
configuration and suppress its speedup.

Speedup uses the fastest validated CPU backend and thread count actually
measured on the same host. A separate speedup relative to the fastest validated
`cpu-batch` configuration is retained.
The scalar tracer lacks separate numerical-failure counters: it enters the
baseline only when it reports no losses and its observables agree with a
no-loss `cpu-batch` run with explicit terminal counts. Collisionless runs
compare moment curves, uncertainty, diffusion and loss count; option 2 also
compares particle histories and transport metadata. These checks establish
agreement on this fixture, not physical convergence. Repeat on each GPU host;
do not divide timings from different hosts. `--compiler` and repeated
`--source-revision NAME=REV` record the build provenance; hardware metadata
and the binary digest are collected automatically. Plotting additionally
requires Matplotlib.

## CINECA Leonardo

The following recipe builds GORILLA_APPLETS and runs the tests on the Leonardo
CPU nodes with GCC 12.2.0, CMake 3.27.9, OpenBLAS 0.3.26, NetCDF-C 4.9.2,
NetCDF-Fortran 4.6.1 and Python 3.11.7. The Spack NetCDF modules ship no
`netCDFConfig.cmake`, and the CMake files link `-lblas` and `-llapack`, which
OpenBLAS does not provide as separate libraries; both are replaced by shims.
`f90nml` for the tests is installed into a virtual environment.

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

Run `ctest` in a batch job on the serial partition
(`--partition=lrd_all_serial --cpus-per-task=4`) rather than on a login node;
the application tests took up to 16 s each in that configuration.
