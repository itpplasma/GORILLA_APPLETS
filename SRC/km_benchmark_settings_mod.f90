module km_benchmark_settings_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    integer, parameter :: max_surfaces = 30
    integer, parameter :: max_background_species = 10
    integer, parameter :: max_trace_scan_points = 12

    integer, public, protected :: collision_operator = 4
    integer, public, protected :: i_integrator_type = 1
    integer, public, protected :: n_particles = 5000
    integer, public, protected :: n_surfaces = 0
    integer, public, protected :: n_background_species = 2
    integer, public, protected :: energy_scan_points = 1
    integer, public, protected :: n_trace_scan_multipliers = 0
    integer, public, protected :: tracer_species = 1
    logical, public, protected :: boole_run_energy_scan = .false.
    logical, public, protected :: boole_run_trace_scan = .false.
    logical, public, protected :: boole_precalc_collisions = .false.
    logical, public, protected :: boole_write_surface_trace = .true.
    real(dp), public, protected :: temperature_eV = 5.0d2
    real(dp), public, protected :: total_time = 0.0_dp
    real(dp), public, protected :: v_E = 0.0_dp
    real(dp), public, protected :: energy_scan_max_factor = 1.0_dp
    real(dp), public, protected :: energy_scan_min_factor = 1.0_dp
    real(dp), public, protected :: background_density(max_background_species) = 0.0_dp
    real(dp), public, protected :: background_temperature(max_background_species) = 0.0_dp
    real(dp), public, protected :: background_mass_amu(max_background_species) = 0.0_dp
    real(dp), public, protected :: background_charge(max_background_species) = 0.0_dp
    real(dp), public, protected :: fit_end_fraction = 1.0_dp
    real(dp), public, protected :: fit_start_fraction = 0.2_dp
    real(dp), public, protected :: surface_s_values(max_surfaces) = 0.0_dp
    real(dp), public, protected :: trace_scan_multipliers(max_trace_scan_points) = 0.0_dp
    real(dp), public, protected :: trace_time_multiplier = 2.0_dp
    character(len=256), public, protected :: diagnostics_prefix = 'km_trace'
    character(len=256), public, protected :: filename_output = 'km_d11_profile.csv'
    character(len=256), public, protected :: collision_profile_file = ''

    ! Per-surface collision parameters (loaded from profile file)
    integer, public, protected :: n_profile_surfaces = 0
    real(dp), public, protected :: profile_s(max_surfaces) = 0.0_dp
    real(dp), public, protected :: profile_temperature(max_surfaces, max_background_species) = 0.0_dp
    real(dp), public, protected :: profile_tracer_energy_eV(max_surfaces) = 0.0_dp
    real(dp), public, protected :: profile_energy_eV(max_surfaces) = 0.0_dp
    real(dp), public, protected :: profile_density(max_surfaces, max_background_species) = 0.0_dp

    public :: load_km_benchmark_inp
    public :: find_nearest_profile_surface

contains

subroutine load_km_benchmark_inp()

    integer :: inp_unit

    namelist /km_benchmark_nml/ collision_operator, i_integrator_type, &
        n_particles, n_surfaces, n_background_species, tracer_species, &
        boole_precalc_collisions, temperature_eV, total_time, v_E, &
        boole_run_energy_scan, energy_scan_points, energy_scan_min_factor, &
        energy_scan_max_factor, boole_run_trace_scan, &
        n_trace_scan_multipliers, trace_scan_multipliers, &
        boole_write_surface_trace, trace_time_multiplier, fit_start_fraction, &
        fit_end_fraction, diagnostics_prefix, &
        background_density, background_temperature, &
        background_mass_amu, background_charge, &
        surface_s_values, filename_output, collision_profile_file

    open(newunit=inp_unit, file='km_benchmark.inp', status='old', action='read')
    read(inp_unit, nml=km_benchmark_nml)
    close(inp_unit)

    call validate_km_benchmark_settings()
    if (len_trim(collision_profile_file) > 0) call load_collision_profile()

    print *, 'GORILLA_APPLETS: Loaded input data from km_benchmark.inp'

end subroutine load_km_benchmark_inp

subroutine validate_km_benchmark_settings()

    if (n_background_species < 1 .or. n_background_species > max_background_species) then
        error stop 'km_benchmark_settings_mod: invalid n_background_species'
    end if
    if (n_surfaces < 0 .or. n_surfaces > max_surfaces) then
        error stop 'km_benchmark_settings_mod: invalid n_surfaces'
    end if
    if (energy_scan_points < 1) then
        error stop 'km_benchmark_settings_mod: energy_scan_points must be positive'
    end if
    if (energy_scan_min_factor <= 0.0_dp .or. energy_scan_max_factor <= 0.0_dp) then
        error stop 'km_benchmark_settings_mod: energy scan factors must be positive'
    end if
    if (energy_scan_max_factor < energy_scan_min_factor) then
        error stop 'km_benchmark_settings_mod: energy_scan_max_factor < energy_scan_min_factor'
    end if
    if (trace_time_multiplier <= 0.0_dp) then
        error stop 'km_benchmark_settings_mod: trace_time_multiplier must be positive'
    end if
    if (fit_start_fraction < 0.0_dp .or. fit_start_fraction >= 1.0_dp) then
        error stop 'km_benchmark_settings_mod: invalid fit_start_fraction'
    end if
    if (fit_end_fraction <= fit_start_fraction .or. fit_end_fraction > 1.0_dp) then
        error stop 'km_benchmark_settings_mod: invalid fit_end_fraction'
    end if
    if (n_trace_scan_multipliers < 0 .or. n_trace_scan_multipliers > max_trace_scan_points) then
        error stop 'km_benchmark_settings_mod: invalid n_trace_scan_multipliers'
    end if
    if (boole_run_trace_scan .and. n_trace_scan_multipliers < 1) then
        error stop 'km_benchmark_settings_mod: trace scan enabled without multipliers'
    end if

end subroutine validate_km_benchmark_settings

subroutine load_collision_profile()

    integer :: profile_unit, i, n_bg, n_cols, n_s
    character(len=256) :: line
    real(dp) :: values(2 * max_background_species + 2)

    open(newunit=profile_unit, file=trim(collision_profile_file), &
        status='old', action='read')

    ! Skip comment lines starting with #
    do
        read(profile_unit, '(A)') line
        if (line(1:1) /= '#') exit
    end do

    ! First non-comment line: n_surfaces n_background_species
    read(line, *) n_s, n_bg
    if (n_s > max_surfaces .or. n_bg > max_background_species) then
        error stop 'km_benchmark_settings_mod: collision profile exceeds configured limits'
    end if
    n_profile_surfaces = n_s

    ! Skip header line
    read(profile_unit, '(A)') line

    ! Read data:
    ! legacy format: s, n_1..n_N, E_eV
    ! current format: s, n_1..n_N, T_1..T_N, E_ref_eV
    do i = 1, n_s
        read(profile_unit, '(A)') line
        values = 0.0_dp
        n_cols = count_tokens(line)
        read(line, *) values(1:n_cols)
        profile_s(i) = values(1)
        profile_density(i, 1:n_bg) = values(2:1 + n_bg)
        if (n_cols == n_bg + 2) then
            profile_temperature(i, 1:n_bg) = values(n_bg + 2)
            profile_tracer_energy_eV(i) = values(n_bg + 2)
        else if (n_cols == 2 * n_bg + 2) then
            profile_temperature(i, 1:n_bg) = values(2 + n_bg:1 + 2 * n_bg)
            profile_tracer_energy_eV(i) = values(2 * n_bg + 2)
        else
            error stop 'km_benchmark_settings_mod: unsupported collision profile format'
        end if
        profile_energy_eV(i) = profile_tracer_energy_eV(i)
    end do
    close(profile_unit)

    print '(A,I0,A,I0,A)', '  Loaded collision profile: ', &
        n_s, ' surfaces, ', n_bg, ' species'

end subroutine load_collision_profile

integer function count_tokens(line) result(n_tokens)

    character(len=*), intent(in) :: line

    integer :: i
    logical :: in_token

    n_tokens = 0
    in_token = .false.
    do i = 1, len_trim(line)
        if (line(i:i) /= ' ' .and. line(i:i) /= char(9)) then
            if (.not. in_token) n_tokens = n_tokens + 1
            in_token = .true.
        else
            in_token = .false.
        end if
    end do

end function count_tokens

integer function find_nearest_profile_surface(target_s) result(profile_index)

    real(dp), intent(in) :: target_s

    real(dp) :: best_distance, distance
    integer :: i

    if (n_profile_surfaces < 1) error stop 'find_nearest_profile_surface: no profile loaded'

    profile_index = 1
    best_distance = abs(profile_s(1) - target_s)
    do i = 2, n_profile_surfaces
        distance = abs(profile_s(i) - target_s)
        if (distance < best_distance) then
            best_distance = distance
            profile_index = i
        end if
    end do

end function find_nearest_profile_surface

end module km_benchmark_settings_mod
