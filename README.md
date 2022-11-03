# GORILLA
**G**uiding-center **OR**bit **I**ntegration with **L**ocal **L**inearization **A**pproach

![alt text](https://github.com/itpplasma/GORILLA/blob/main/DOCUMENTATION/LATEX/figures/title_image.png "Guiding-center orbit Poincaré plot with GORILLA")

GORILLA computes guiding-center orbits for charged particles of given mass, charge and energy in toroidal fusion devices with three-dimensional field geometry. This guiding-center orbit code is designed to be used in scientific plasma physics simulations in the field of magnetic confinement fusion.

### Summary
The guiding-center orbits are traced via a quasi-geometric integration method described in Ref. [1].
There, high order interpolation of electromagnetic fields in space is replaced by a special linear interpolation, leading to locally linear Hamiltonian equations of motion with piecewise constant coefficients. The underlying formulation treats the motion in the piecewise linear fields exactly. This further leads to conservation of total energy, magnetic moment and phase space volume. Furthermore, the approach reduces computational effort and noise sensitivity. Guiding-center orbits are computed without taking collisions into account.

Due to its formulation in general curvilinear coordinates, GORILLA is not limited by the field topology. That means that the computation domain of GORILLA covers both the closed field line region (i.e. the plasma core) and the open field line region (i.e. the scrape-off layer).

For various simulations in magnetic confinement fusion, direct modeling of guiding-center particle orbits is utilized, e.g. global kinetic computations of quasi-steady plasma parameters or fast alpha particle loss estimation for stellarator optimization. In such complex simulations a simple interface for the guiding-center orbit integration part is needed. Namely, the initial condition in five-dimensional phase space is provided (i.e. guiding-center position, parallel and perpendicular velocity) and the main interest is in the condition after a prescribed time step while the integration process itself is irrelevant. Such a pure “orbit time step routine” acting as an interface with a plasma physics simulation is provided (`orbit_timestep_gorilla`).
However, the integration process itself can be of high interest as well, thus, a program allowing the detailed analysis of guiding-center orbits, the time evolution of their respective invariants of motion and Poincaré plots is at disposal as well (`gorilla_plot`).
Both applications are realized for demonstration in the program (`test_gorilla_main`).

### License
GORILLA_APPLETS is mainly intended to be used privatly by the Theoretical Plasma group of the Theoretical and Computational Physics Departement of the Technical University Graz.


## Documentation
The following supplemental material is available in `DOCUMENTATION/SUPPLEMENTAL_MATERIAL`:
* Master's thesis of J. Schatzlmayr (Bolzmann test, dwell time calculation)
* Master's thesis of G. Graßler (reversibility test)


### Building: How to work with the main GORILLA


### Building with make
```bash
cd /path/to/GORILLA
make
```
This will produce `gorilla_applets_main.x` required to run the code. To specify the location of
NetCDF includes and libraries, one has to set the `NCINC` and `NCLIB` variable during `make`.


## Usage

GORILLA and therefore GORILLA_APPLETS currently runs on a single node with OpenMP shared memory parallelization with one particle per thread and background fields residing in main memory.

The main executable is `gorilla_applets_main.x`.
As an input it takes ....

### ... of the main GORILLA (see respective folder/repository)
... the following input files
* `tetra_grid.inp`                       (Input file for settings of the tetrahedronal grid used in GORLLA)
* `gorilla.inp`                             (Input file for settings of GORILLA)
* `gorilla_plot.inp`                   (Input file for the program for the analysis of guiding-center orbits)
* `field_divB0.inp`                     (Input file for loading g-file equilibria - Do not change this file.)
* `preload_for_SYNCH.inp`         (Input file for splining magnetic field data of g-file equilibria - Do not change this file.)

... and the MHD equilibrium files 
* `netcdf_file_for_test.nc`: VMEC NetCDF equlibrium (File name can be specified in `tetra_grid.inp`.)
* `g_file_for_test` or `g_file_for_test_WEST`: g-file equilibrium (File name can be specified in `tetra_grid.inp`.)

For compability with WEST geometry of SOLEDGE3X-EIRENE, additional input files describing the original 2D mesh are needed. 

* `knots_for_test.dat`: coordinates ($R$, $Z$) of the vertices making up the 2D grid (File name can be specified in `tetra_grid.inp`.)
* `triangles_for_test.dat`: association of above mentioned vertices to triangles (triples of vertices) covering the 2D plane (File name can be specified in `tetra_grid.inp`.)

### ... specifically for GORILLA_APPLETS
... the following input files which can be found in the folder `INPUT/`
* `gorilla_applets.inp`                       (Input file for main setting of GORILLA_APPLETS like choosing the application and precomputation of fluxtubevolumns)
* `mono_energetic_tranps_coef.inp`                             (Input file for settings of calculation of monoenergetic transport coefficient)
* `alpha_lifetime.inp`                   (Input file for the program for the analysis of the alpha particle loss)
* `total_dwell_times.inp`                     (Input file for the calculation of particle dwell times in the tetrahedral cells.)
* `boltzmann.inp`         (Input file for the bolzmann test.)
* `reversibility_test.inp`         (Input file for the reversibility test, demonstrating artifact-free integration.)
* `seed.inp`         (Example seed used for several instances of random number generation - Can be generated by GORILLA, see respective .inp file of the application in question.)

## Examples

Three examples for calculating alpha particle loss, demonstrating the conservation of the poincaré invariant and showing artifact-free integration via an reversibility test can be found in `EXAMPLES/`. There, the necessary soft links are already created and the input files are given, and runs are started with
```
./gorilla_applets_main.x   #if the build was done with make
```
To avoid hyperthreading issues, it is beneficial to limit the number of threads to
the number of actual CPU cores via the environment variable `$OMP_NUM_THREADS`.
Detailed descriptions of the respective input files specific to GORILLA_APPLETS can be found in `INPUT`.
After appropriate compilation of GORILLA_APPLETS, the code can be executed in all of these 3 example folders, respectively.
For the visualization of the output of these three examples, appropriate plotting methods for MATLAB are at disposal in the respective folders.

### ALPHA_LIFETIME
* 

### POINCARE_INVARIANT
* 

### REVERSIBILITY_TEST
* 


### Generation of input files and plotting in MATLAB and Python
A detailed explanation of all examples (1-7) including the generation of the appropriate input files (including the example folders in `EXAMPLES/MATLAB_RUN` and `EXAMPLES/PYTHON_RUN`) and plotting of the results with MATLAB and Python can be found in the folders `MATLAB` and `PYTHON`, respectively.
Here, the results of GORILLA with different polynominal orders K=2,3,4 and Runge-Kutta 4 are compared in case of examples 1-3. For examples 5-6 orbits for both trapped and passing particles are calculated. For example 7 an additional, in-depth comparison between adaptive and non-adaptive scheme is performed.


## References for GORILLA
When using this code for scientific publications, please cite both references [1] and [2]:

[1] M. Eder, C.G. Albert, L.M.P. Bauer, S.V. Kasilov and W. Kernbichler
“Quasi-geometric integration of guiding-center orbits in piecewise linear toroidal fields”
Physics of Plasmas 27, 122508 (2020)
<https://doi.org/10.1063/5.0022117>
Preprint: <https://arxiv.org/abs/2007.08151>

[2] M. Eder, C.G. Albert, L.M.P. Bauer, Georg S. Graßler, S.V. Kasilov, W. Kernbichler, M. Meisterhofer and M.Scheidt
“GORILLA: Guiding-center ORbit Integration with Local Linearization Approach”
submitted to Journal of Open Source Software
Preprint: <https://github.com/openjournals/joss-papers/blob/joss.03116/joss.03116/10.21105.joss.03116.pdf>
