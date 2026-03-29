module usual_transport_benchmark_settings_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    integer, public, protected :: collision_operator = 4
    integer, public, protected :: i_integrator_type = 1
    integer, public, protected :: n_particles = 10000
    integer, public, protected :: tracer_species = 1
    logical, public, protected :: boole_precalc_collisions = .false.
    real(dp), public, protected :: density = 0.0_dp
    real(dp), public, protected :: density_log_gradient_per_s = 0.0_dp
    real(dp), public, protected :: density_profile_reference_s = 0.5_dp
    real(dp), public, protected :: electron_density = 0.0_dp
    real(dp), public, protected :: electron_temperature_eV = 0.0_dp
    real(dp), public, protected :: ion_density = 0.0_dp
    real(dp), public, protected :: ion_temperature_eV = 0.0_dp
    real(dp), public, protected :: surface_s = 0.5_dp
    real(dp), public, protected :: temperature_eV = 0.0_dp
    real(dp), public, protected :: total_time = 1.0e-4_dp
    real(dp), public, protected :: v_E = 0.0_dp
    character(len=128), public, protected :: filename_transport_summary = 'usual_transport_summary.dat'
    character(len=128), public, protected :: filename_boundary_flux = 'usual_transport_boundary_flux.dat'

    public :: load_usual_transport_benchmark_inp

contains

subroutine load_usual_transport_benchmark_inp()

    integer :: benchmark_unit

    namelist /usual_transport_benchmark_nml/ i_integrator_type, n_particles, tracer_species, boole_precalc_collisions, &
        density, density_log_gradient_per_s, density_profile_reference_s, electron_density, &
        electron_temperature_eV, ion_density, ion_temperature_eV, surface_s, temperature_eV, total_time, v_E, &
        collision_operator, filename_transport_summary, filename_boundary_flux

    open(newunit=benchmark_unit, file='usual_transport_benchmark.inp', status='old', action='read')
    read(benchmark_unit, nml=usual_transport_benchmark_nml)
    close(benchmark_unit)

    if (electron_density <= 0.0_dp) electron_density = density
    if (ion_density <= 0.0_dp) ion_density = density
    if (electron_temperature_eV <= 0.0_dp) electron_temperature_eV = temperature_eV
    if (ion_temperature_eV <= 0.0_dp) ion_temperature_eV = temperature_eV
    if (density_profile_reference_s < 0.0_dp) density_profile_reference_s = surface_s

    print *, 'GORILLA_APPLETS: Loaded input data from usual_transport_benchmark.inp'

end subroutine load_usual_transport_benchmark_inp

end module usual_transport_benchmark_settings_mod
