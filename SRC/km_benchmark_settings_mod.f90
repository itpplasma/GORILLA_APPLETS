module km_benchmark_settings_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    integer, parameter :: max_surfaces = 30

    integer, public, protected :: collision_operator = 4
    integer, public, protected :: i_integrator_type = 1
    integer, public, protected :: n_particles = 5000
    integer, public, protected :: n_surfaces = 0
    integer, public, protected :: tracer_species = 1
    logical, public, protected :: boole_precalc_collisions = .false.
    real(dp), public, protected :: electron_density = 3.0d13
    real(dp), public, protected :: electron_temperature_eV = 3.5d3
    real(dp), public, protected :: ion_density = 3.0d13
    real(dp), public, protected :: ion_temperature_eV = 3.5d3
    real(dp), public, protected :: temperature_eV = 3.5d3
    real(dp), public, protected :: total_time = 0.0_dp
    real(dp), public, protected :: v_E = 0.0_dp
    real(dp), public, protected :: surface_s_values(max_surfaces) = 0.0_dp
    character(len=256), public, protected :: filename_output = 'km_d11_profile.csv'

    public :: load_km_benchmark_inp

contains

subroutine load_km_benchmark_inp()

    integer :: inp_unit

    namelist /km_benchmark_nml/ collision_operator, i_integrator_type, n_particles, n_surfaces, &
        tracer_species, boole_precalc_collisions, electron_density, electron_temperature_eV, &
        ion_density, ion_temperature_eV, temperature_eV, total_time, v_E, surface_s_values, &
        filename_output

    open(newunit=inp_unit, file='km_benchmark.inp', status='old', action='read')
    read(inp_unit, nml=km_benchmark_nml)
    close(inp_unit)

    if (electron_density <= 0.0_dp) electron_density = 3.0d13
    if (ion_density <= 0.0_dp) ion_density = 3.0d13
    if (electron_temperature_eV <= 0.0_dp) electron_temperature_eV = temperature_eV
    if (ion_temperature_eV <= 0.0_dp) ion_temperature_eV = temperature_eV

    print *, 'GORILLA_APPLETS: Loaded input data from km_benchmark.inp'

end subroutine load_km_benchmark_inp

end module km_benchmark_settings_mod
