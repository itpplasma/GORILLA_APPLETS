module km_benchmark_settings_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    integer, parameter :: max_surfaces = 30
    integer, parameter :: max_background_species = 10

    integer, public, protected :: collision_operator = 4
    integer, public, protected :: i_integrator_type = 1
    integer, public, protected :: n_particles = 5000
    integer, public, protected :: n_surfaces = 0
    integer, public, protected :: n_background_species = 2
    integer, public, protected :: tracer_species = 1
    logical, public, protected :: boole_precalc_collisions = .false.
    real(dp), public, protected :: temperature_eV = 5.0d2
    real(dp), public, protected :: total_time = 0.0_dp
    real(dp), public, protected :: v_E = 0.0_dp
    real(dp), public, protected :: background_density(max_background_species) = 0.0_dp
    real(dp), public, protected :: background_temperature(max_background_species) = 0.0_dp
    real(dp), public, protected :: background_mass_amu(max_background_species) = 0.0_dp
    real(dp), public, protected :: background_charge(max_background_species) = 0.0_dp
    real(dp), public, protected :: surface_s_values(max_surfaces) = 0.0_dp
    character(len=256), public, protected :: filename_output = 'km_d11_profile.csv'

    public :: load_km_benchmark_inp

contains

subroutine load_km_benchmark_inp()

    integer :: inp_unit

    namelist /km_benchmark_nml/ collision_operator, i_integrator_type, &
        n_particles, n_surfaces, n_background_species, tracer_species, &
        boole_precalc_collisions, temperature_eV, total_time, v_E, &
        background_density, background_temperature, &
        background_mass_amu, background_charge, &
        surface_s_values, filename_output

    open(newunit=inp_unit, file='km_benchmark.inp', status='old', action='read')
    read(inp_unit, nml=km_benchmark_nml)
    close(inp_unit)

    print *, 'GORILLA_APPLETS: Loaded input data from km_benchmark.inp'

end subroutine load_km_benchmark_inp

end module km_benchmark_settings_mod
