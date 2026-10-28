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

The executable is `BUILD/gorilla_applets_main.x`. The example folders expect
a link to it in the repository root (`ln -s BUILD/gorilla_applets_main.x .`).

## Tests

The CTest suite consists of Fortran unit tests (collision operators, RMP
loading and profile moments, transport normalisation and loss moments) and
short application runs driven by the Python scripts in `TESTS/`, for example
the alpha-particle lifetime test, which runs the flux-tube volume
precomputation and then a ten-particle trace and checks the output.
`TESTS/OU_COLLISION_INVARIANT.md` and `TESTS/RMP_DENSITY_RESPONSE.md` describe
two of the physics checks.

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
