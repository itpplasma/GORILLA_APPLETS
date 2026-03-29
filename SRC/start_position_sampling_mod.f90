module start_position_sampling_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    public :: set_random_start_positions
    public :: validate_current_start_positions

contains

subroutine set_random_start_positions(rand_matrix, species_in, s0, s_min_in, s_max_in)

    use gorilla_applets_types_mod, only: in, start
    use tetra_grid_settings_mod, only: sfc_s_min

    real(dp), intent(in) :: rand_matrix(:, :, :)
    integer, intent(in), optional :: species_in(:)
    real(dp), intent(in), optional :: s0
    real(dp), intent(in), optional :: s_min_in
    real(dp), intent(in), optional :: s_max_in

    integer, allocatable :: species(:)
    integer :: i
    integer :: primary_species
    logical :: use_fixed_s
    real(dp) :: fixed_s
    real(dp) :: s_max
    real(dp) :: s_min

    if (present(species_in)) then
        allocate(species(size(species_in)))
        species = species_in
    else
        allocate(species(in%n_species))
        do i = 1, in%n_species
            species(i) = i
        end do
    end if

    s_min = sfc_s_min
    s_max = 1.0_dp
    if (present(s_min_in)) s_min = max(sfc_s_min, min(s_min_in, 1.0_dp))
    if (present(s_max_in)) s_max = max(s_min, min(s_max_in, 1.0_dp))

    use_fixed_s = present(s0)
    if (use_fixed_s) fixed_s = s0

    primary_species = species(1)
    call fill_random_start_positions(rand_matrix, primary_species, s_min, s_max, use_fixed_s, fixed_s)

    if (size(species) > 1) then
        do i = 2, size(species)
            start%x(:, :, species(i)) = start%x(:, :, primary_species)
        end do
    end if

end subroutine set_random_start_positions

subroutine validate_current_start_positions(species_in, s0, s_min_in, s_max_in, max_validation_attempts, failed_positions)

    use gorilla_applets_types_mod, only: in, start
    use tetra_grid_settings_mod, only: sfc_s_min

    integer, intent(in), optional :: species_in(:)
    real(dp), intent(in), optional :: s0
    real(dp), intent(in), optional :: s_min_in
    real(dp), intent(in), optional :: s_max_in
    integer, intent(in), optional :: max_validation_attempts
    integer, intent(out), optional :: failed_positions

    integer, allocatable :: species(:)
    integer :: failed_count
    integer :: i
    integer :: max_attempts
    integer :: primary_species
    logical :: use_fixed_s
    real(dp) :: fixed_s
    real(dp) :: s_max
    real(dp) :: s_min

    if (present(species_in)) then
        allocate(species(size(species_in)))
        species = species_in
    else
        allocate(species(in%n_species))
        do i = 1, in%n_species
            species(i) = i
        end do
    end if

    s_min = sfc_s_min
    s_max = 1.0_dp
    if (present(s_min_in)) s_min = max(sfc_s_min, min(s_min_in, 1.0_dp))
    if (present(s_max_in)) s_max = max(s_min, min(s_max_in, 1.0_dp))

    use_fixed_s = present(s0)
    if (use_fixed_s) fixed_s = s0

    max_attempts = 64
    if (present(max_validation_attempts)) max_attempts = max(1, max_validation_attempts)

    primary_species = species(1)
    call validate_start_positions(primary_species, s_min, s_max, use_fixed_s, fixed_s, max_attempts, failed_count)

    if (size(species) > 1) then
        do i = 2, size(species)
            start%x(:, :, species(i)) = start%x(:, :, primary_species)
        end do
    end if

    if (present(failed_positions)) failed_positions = failed_count

end subroutine validate_current_start_positions

subroutine fill_random_start_positions(rand_matrix, species, s_min, s_max, use_fixed_s, fixed_s)

    use constants, only: pi
    use gorilla_applets_types_mod, only: start
    use tetra_grid_settings_mod, only: n_field_periods

    real(dp), intent(in) :: rand_matrix(:, :, :)
    integer, intent(in) :: species
    logical, intent(in) :: use_fixed_s
    real(dp), intent(in) :: fixed_s
    real(dp), intent(in) :: s_max
    real(dp), intent(in) :: s_min

    start%x(1, :, species) = s_min + rand_matrix(1, :, 1) * (s_max - s_min)
    if (use_fixed_s) start%x(1, :, species) = fixed_s
    start%x(2, :, species) = 2.0_dp * pi * rand_matrix(2, :, 1)
    start%x(3, :, species) = 2.0_dp * pi * rand_matrix(3, :, 1) / real(n_field_periods, dp)

end subroutine fill_random_start_positions

subroutine validate_start_positions(species, s_min, s_max, use_fixed_s, fixed_s, max_attempts, failed_positions)

    use constants, only: ev2erg, pi
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_types_mod, only: start
    use tetra_grid_settings_mod, only: n_field_periods

    integer, intent(in) :: species
    logical, intent(in) :: use_fixed_s
    integer, intent(in) :: max_attempts
    integer, intent(out) :: failed_positions
    real(dp), intent(in) :: fixed_s
    real(dp), intent(in) :: s_max
    real(dp), intent(in) :: s_min

    integer :: attempt
    integer :: iface
    integer :: ind_tetr
    integer :: n
    real(dp) :: random_values(3)
    real(dp) :: speed
    real(dp) :: vpar
    real(dp) :: vperp
    real(dp) :: x(3)

    failed_positions = 0

    do n = 1, size(start%x, 2)
        x = start%x(:, n, species)
        speed = sqrt(2.0_dp * start%energy(n, species) * ev2erg / start%particle_mass(species))
        vpar = start%pitch(n, species) * speed
        vperp = sqrt(max(speed * speed - vpar * vpar, 0.0_dp))
        call find_tetra(x, vpar, vperp, ind_tetr, iface)
        attempt = 0
        do while (ind_tetr == -1 .and. attempt < max_attempts)
            attempt = attempt + 1
            call random_number(random_values)
            x(1) = s_min + random_values(1) * (s_max - s_min)
            if (use_fixed_s) x(1) = fixed_s
            x(2) = 2.0_dp * pi * random_values(2)
            x(3) = 2.0_dp * pi * random_values(3) / real(n_field_periods, dp)
            call find_tetra(x, vpar, vperp, ind_tetr, iface)
        end do
        if (ind_tetr == -1) failed_positions = failed_positions + 1
        start%x(:, n, species) = x
    end do

end subroutine validate_start_positions

end module start_position_sampling_mod
