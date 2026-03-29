module global_transport_fit_io_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    use global_transport_fit_math_mod, only: compute_geometry_from_boundaries
    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_fit_result_t

    implicit none

    private

    public :: load_global_transport_experiments
    public :: write_convergence_history
    public :: write_global_transport_fit_outputs

contains

subroutine load_global_transport_experiments(boundary_file, shell_volume_file, source_files, density_files, density_var_files, &
    flux_files, flux_var_files, experiments)

    character(len=*), intent(in) :: boundary_file
    character(len=*), intent(in) :: shell_volume_file
    character(len=*), intent(in) :: source_files(:)
    character(len=*), intent(in) :: density_files(:)
    character(len=*), intent(in) :: density_var_files(:)
    character(len=*), intent(in) :: flux_files(:)
    character(len=*), intent(in) :: flux_var_files(:)
    type(global_transport_experiment_t), allocatable, intent(out) :: experiments(:)

    real(dp), allocatable :: boundary_s(:)
    real(dp), allocatable :: boundary_areas(:)
    real(dp), allocatable :: cell_centers(:)
    real(dp), allocatable :: shell_volumes(:)
    integer :: i

    call read_vector_file(boundary_file, boundary_s)
    call read_vector_file(shell_volume_file, shell_volumes)
    call compute_geometry_from_boundaries(boundary_s, shell_volumes, cell_centers, boundary_areas)
    allocate(experiments(size(source_files)))

    do i = 1, size(source_files)
        experiments(i)%boundary_s = boundary_s
        experiments(i)%shell_volumes = shell_volumes
        call read_vector_file(source_files(i), experiments(i)%source)
        call read_vector_file(density_files(i), experiments(i)%density)
        call read_vector_file(density_var_files(i), experiments(i)%density_variance)
        call read_vector_file(flux_files(i), experiments(i)%flux)
        call read_vector_file(flux_var_files(i), experiments(i)%flux_variance)
        experiments(i)%integrated_flux = boundary_areas * experiments(i)%flux
        experiments(i)%integrated_flux_variance = boundary_areas**2 * experiments(i)%flux_variance
    end do

end subroutine load_global_transport_experiments

subroutine write_global_transport_fit_outputs(boundary_s, result, summary_file, profiles_file)

    real(dp), intent(in) :: boundary_s(:)
    type(global_transport_fit_result_t), intent(in) :: result
    character(len=*), intent(in) :: summary_file
    character(len=*), intent(in) :: profiles_file

    integer :: io_unit
    integer :: i

    open(newunit=io_unit, file=summary_file, status='replace', action='write')
    write(io_unit, '(A)') 'objective converged iterations n_parameters'
    write(io_unit, '(ES24.16,1X,L1,1X,I0,1X,I0)') result%objective, result%converged, result%n_iterations, size(result%parameters)
    write(io_unit, '(A)') 'parameter_index parameter gradient variance'
    do i = 1, size(result%parameters)
        write(io_unit, '(I0,3(1X,ES24.16))') i, result%parameters(i), result%gradient(i), result%covariance(i, i)
    end do
    close(io_unit)

    open(newunit=io_unit, file=profiles_file, status='replace', action='write')
    write(io_unit, '(A)') 'boundary_index boundary_s A_s A_s_std A_s_2sigma B_s B_s_std B_s_2sigma'
    do i = 1, size(boundary_s)
        write(io_unit, '(I0,7(1X,ES24.16))') i, boundary_s(i), result%a_profile(i), result%a_std(i), &
            2.0_dp * result%a_std(i), result%b_profile(i), result%b_std(i), 2.0_dp * result%b_std(i)
    end do
    close(io_unit)

end subroutine write_global_transport_fit_outputs

subroutine write_convergence_history(result, filename)

    type(global_transport_fit_result_t), intent(in) :: result
    character(len=*), intent(in) :: filename

    integer :: i
    integer :: io_unit

    if (.not.allocated(result%history_objective)) return

    open(newunit=io_unit, file=filename, status='replace', action='write')
    write(io_unit, '(A)') 'iteration objective gradient_norm step_norm damping accepted'
    do i = 1, size(result%history_objective)
        write(io_unit, '(I0,4(1X,ES24.16),1X,L1)') i, result%history_objective(i), &
            result%history_gradient_norm(i), result%history_step_norm(i), &
            result%history_damping(i), result%history_accepted(i)
    end do
    close(io_unit)

end subroutine write_convergence_history

subroutine read_vector_file(filename, vector)

    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: vector(:)

    integer :: count
    integer :: io_unit
    integer :: io_status
    real(dp) :: value
    integer :: i

    count = 0
    open(newunit=io_unit, file=filename, status='old', action='read')
    do
        read(io_unit, *, iostat=io_status) value
        if (io_status /= 0) exit
        count = count + 1
    end do
    rewind(io_unit)

    allocate(vector(count))
    do i = 1, count
        read(io_unit, *) vector(i)
    end do
    close(io_unit)

end subroutine read_vector_file

end module global_transport_fit_io_mod
