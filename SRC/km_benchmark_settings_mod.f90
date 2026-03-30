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
    character(len=256), public, protected :: collision_profile_file = ''

    ! Per-surface collision parameters (loaded from profile file)
    integer, public, protected :: n_profile_surfaces = 0
    real(dp), public, protected :: profile_s(max_surfaces) = 0.0_dp
    real(dp), public, protected :: profile_energy_eV(max_surfaces) = 0.0_dp
    real(dp), public, protected :: profile_density(max_surfaces, max_background_species) = 0.0_dp

    public :: load_km_benchmark_inp

contains

subroutine load_km_benchmark_inp()

    integer :: inp_unit

    namelist /km_benchmark_nml/ collision_operator, i_integrator_type, &
        n_particles, n_surfaces, n_background_species, tracer_species, &
        boole_precalc_collisions, temperature_eV, total_time, v_E, &
        background_density, background_temperature, &
        background_mass_amu, background_charge, &
        surface_s_values, filename_output, collision_profile_file

    open(newunit=inp_unit, file='km_benchmark.inp', status='old', action='read')
    read(inp_unit, nml=km_benchmark_nml)
    close(inp_unit)

    if (len_trim(collision_profile_file) > 0) call load_collision_profile()

    print *, 'GORILLA_APPLETS: Loaded input data from km_benchmark.inp'

end subroutine load_km_benchmark_inp

subroutine load_collision_profile()

    integer :: profile_unit, i, n_s, n_bg
    character(len=256) :: line

    open(newunit=profile_unit, file=trim(collision_profile_file), &
        status='old', action='read')

    ! Skip comment lines starting with #
    do
        read(profile_unit, '(A)') line
        if (line(1:1) /= '#') exit
    end do

    ! First non-comment line: n_surfaces n_background_species
    read(line, *) n_s, n_bg
    n_profile_surfaces = n_s

    ! Skip header line
    read(profile_unit, '(A)') line

    ! Read data: s, n_1..n_N, E_eV
    do i = 1, n_s
        read(profile_unit, *) profile_s(i), &
            profile_density(i, 1:n_bg), profile_energy_eV(i)
    end do
    close(profile_unit)

    print '(A,I0,A,I0,A)', '  Loaded collision profile: ', &
        n_s, ' surfaces, ', n_bg, ' species'

end subroutine load_collision_profile

end module km_benchmark_settings_mod
