module boltzmann_writing_data_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    public :: write_data_to_files, name_files, unlink_files

contains

subroutine write_data_to_files(filenames,output,moment_specs)

    use boltzmann_types_mod, only: filenames_t, boole_writing_data_t, output_t, moment_specs_t

    type(filenames_t), intent(in) :: filenames
    type(output_t), intent(in) :: output
    type(moment_specs_t), intent(in) :: moment_specs
    type(boole_writing_data_t) :: boole_writing_data

    if (boole_writing_data%vertex_indices) &
        call write_vertex_indices(filenames%vertex_indices)

    if (boole_writing_data%vertex_coordinates) &
        call write_vertex_coordinates(filenames%vertex_coordinates)

    if (boole_writing_data%prism_volumes) &
        call write_prism_volumes(filenames%prism_volumes,output%prism_volumes)

    if (boole_writing_data%refined_prism_volumes) &
        call write_refined_prism_volumes(filenames%refined_prism_volumes,output%refined_prism_volumes)

    if (boole_writing_data%boltzmann_density) &
        call write_boltzmann_densities(filenames%boltzmann_density,output%boltzmann_density)

    if (boole_writing_data%electric_potential) &
        call write_electric_potential(filenames%electric_potential,output%electric_potential)

    if (boole_writing_data%moments) then
        if (moment_specs%n_moments.gt.0) then
            call write_moments(filenames%prism_moments,filenames%prism_moments_summed_squares,filenames%tetr_moments,moment_specs, &
                               output%prism_moments,output%prism_moments_squared,output%tetr_moments)
        else
            print*, "Error: moments are not written to file because no moment was computed. Turn computation of moments on in &
                     gorilla.inp."
        endif
    endif

    if (boole_writing_data%fourier_moments) then
        if (moment_specs%n_moments.gt.0) then
            call write_fourier_moments(filenames%fourier_moments,moment_specs, output%moments_in_frequency_space)
        else
            print*, "Error: Fourier moments are not written to file because no moment was computed (and thus also no fourier &
                     moment). Turn computation of moments on in gorilla.inp."
        endif

    endif

end subroutine write_data_to_files

subroutine write_vertex_indices(filename_vertex_indices)

    use tetra_grid_mod, only: ntetr, tetra_grid

    character(len=100) :: filename_vertex_indices
    integer :: vi_unit
    integer :: i

    open(newunit = vi_unit, file = filename_vertex_indices)
    do i=1, ntetr
        write(vi_unit, *) tetra_grid(i)%ind_knot([1, 2, 3, 4])
    end do
    close(vi_unit)

end subroutine write_vertex_indices

subroutine write_vertex_coordinates(filename_vertex_coordinates)

    use tetra_physics_mod, only: coord_system
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, nvert

    character(len=100) :: filename_vertex_coordinates
    integer :: vc_unit
    integer :: i

    open(newunit = vc_unit, file = filename_vertex_coordinates)
    101 format(1000(e21.14,x))
    if (coord_system.eq.1) then
        ![R,phi,Z]: Write vertex coordinates to file
        do i=1, nvert
            write(vc_unit,101) verts_rphiz(1, i), verts_rphiz(2, i), verts_rphiz(3, i)
        end do
    elseif (coord_system.eq.2) then
        ![s,theta,phi]: Write vertex coordinates to file
        do i=1, nvert
            write(vc_unit,101) verts_sthetaphi(1, i), verts_sthetaphi(2, i), verts_sthetaphi(3, i)
        end do
    endif
    close(vc_unit)

end subroutine write_vertex_coordinates

subroutine write_prism_volumes(filename_prism_volumes,prism_volumes)

    real(dp),dimension(:), intent(in) :: prism_volumes
    character(len=100) :: filename_prism_volumes
    integer :: pv_unit

    open(newunit = pv_unit, file = filename_prism_volumes)
    write(pv_unit,'(ES20.10E4)') prism_volumes
    close(pv_unit)

end subroutine write_prism_volumes

subroutine write_refined_prism_volumes(filename_refined_prism_volumes,refined_prism_volumes)

    real(dp),dimension(:), intent(in) :: refined_prism_volumes
    character(len=100) :: filename_refined_prism_volumes
    integer :: rpv_unit

    open(newunit = rpv_unit, file = filename_refined_prism_volumes)
    write(rpv_unit,'(ES20.10E4)') refined_prism_volumes
    close(rpv_unit)

end subroutine write_refined_prism_volumes

subroutine write_boltzmann_densities(filename_boltzmann_density,boltzmann_density)

    real(dp),dimension(:), intent(in) :: boltzmann_density
    character(len=100) :: filename_boltzmann_density
    integer :: bd_unit

    open(newunit = bd_unit, file = filename_boltzmann_density)
    write(bd_unit,'(ES20.10E4)') boltzmann_density
    close(bd_unit)

end subroutine write_boltzmann_densities

subroutine write_electric_potential(filename_electric_potential,electric_potential)

    real(dp),dimension(:), intent(in) :: electric_potential
    character(len=100) :: filename_electric_potential
    integer :: epv_unit

    open(newunit = epv_unit, file = filename_electric_potential)
    write(epv_unit,'(ES20.10E4)') electric_potential
    close(epv_unit)

end subroutine write_electric_potential

subroutine write_moments(filename_prism_moments, filename_prism_moments_summed_squares, filename_tetr_moments, moment_specs, &
                         prism_moments, prism_moments_squared, tetr_moments)

    use tetra_grid_mod, only: ntetr
    use boltzmann_types_mod, only: moment_specs_t

    complex(dp), dimension(:,:), intent(in) :: tetr_moments, prism_moments, prism_moments_squared
    type(moment_specs_t), intent(in) :: moment_specs
    character(len=100) :: filename_prism_moments, filename_prism_moments_summed_squares, filename_tetr_moments
    integer :: p_moments_unit, pmss_unit, t_moments_unit
    integer :: l, i
    integer :: n_prisms

    open(newunit = p_moments_unit, file = filename_prism_moments)
    open(newunit = pmss_unit, file = filename_prism_moments_summed_squares)
    open(newunit = t_moments_unit, file = filename_tetr_moments)

    n_prisms = ntetr/3

    if (moment_specs%n_moments.gt.0) then
        do l = 1,n_prisms
            do i = 1,moment_specs%n_moments - 1
                write(p_moments_unit,'(2ES20.10E4)',advance="no") real(prism_moments(i,l)), aimag(prism_moments(i,l))
            enddo
                write(p_moments_unit,'(2ES20.10E4)') real(prism_moments(moment_specs%n_moments,l)), &
                                                     aimag(prism_moments(moment_specs%n_moments,l))
        enddo
        if (moment_specs%boole_squared_moments) then
            do l = 1,n_prisms
                do i = 1,moment_specs%n_moments - 1
                    write(pmss_unit,'(2ES20.10E4)',advance="no") real(prism_moments_squared(i,l)), &
                                                                    & aimag(prism_moments_squared(i,l))
                enddo
                    write(pmss_unit,'(2ES20.10E4)') real(prism_moments_squared(moment_specs%n_moments,l)), &
                                                    aimag(prism_moments_squared(moment_specs%n_moments,l))
            enddo
        endif
        do l = 1,ntetr
            do i = 1,moment_specs%n_moments - 1
                write(t_moments_unit,'(2ES20.10E4)',advance="no") real(tetr_moments(i,l)), aimag(tetr_moments(i,l))
            enddo
                write(t_moments_unit,'(2ES20.10E4)') real(tetr_moments(moment_specs%n_moments,l)), &
                                                     aimag(tetr_moments(moment_specs%n_moments,l))
        enddo
    endif

    close(p_moments_unit)
    close(pmss_unit)
    close(t_moments_unit)

end subroutine write_moments

subroutine write_fourier_moments(filename_fourier_moments,moment_specs,moments_in_frequency_space)

    use tetra_grid_settings_mod, only: grid_size
    use boltzmann_types_mod, only: moment_specs_t

    complex(dp), dimension(:,:,:), intent(in) :: moments_in_frequency_space
    type(moment_specs_t), intent(in) :: moment_specs
    character(len=100) :: filename_fourier_moments
    integer :: fm_unit, l, i

    open(newunit = fm_unit, file = filename_fourier_moments)
    do l = 1,moment_specs%n_triangles
        do i = 1,moment_specs%n_fourier_modes-1
            write(fm_unit,'(2ES20.10E4)',advance="no") real(moments_in_frequency_space(1,l,i)), &
                                                       aimag(moments_in_frequency_space(1,l,i))
        enddo
            write(fm_unit,'(2ES20.10E4)') real(moments_in_frequency_space(1,l,moment_specs%n_fourier_modes)), &
                                          aimag(moments_in_frequency_space(1,l,moment_specs%n_fourier_modes))
    enddo
    close(fm_unit)

end subroutine write_fourier_moments

subroutine name_files(filenames)

    use boltzmann_types_mod, only: filenames_t

    type(filenames_t), intent(out) :: filenames

    filenames%exit_times = 'exit_times.dat'
    filenames%remaining_particles = 'remaining_particles.dat'
    filenames%poincare_maps = 'poincare_maps.dat'
    filenames%prism_moments = 'prism_moments.dat'
    filenames%prism_moments_summed_squares = 'prism_moments_summed_squares.dat'
    filenames%vertex_coordinates = 'vertex_coordinates.dat'
    filenames%vertex_indices = 'vertex_indices.dat'
    filenames%prism_volumes = 'prism_volumes.dat'
    filenames%fourier_moments = 'fourier_moments.dat'
    filenames%refined_prism_volumes = 'refined_prism_volumes.dat'
    filenames%electric_potential = 'electric_potential.dat'
    filenames%boltzmann_density = 'boltzmann_density.dat'
    filenames%divertor_intersections = 'divertor_intersections.dat'
    filenames%tetr_moments = 'tetr_moments.dat'

end subroutine name_files

subroutine unlink_files(filenames)

    use boltzmann_types_mod, only: filenames_t

    type(filenames_t) :: filenames

    call unlink(filenames%exit_times)
    call unlink(filenames%remaining_particles)
    call unlink(filenames%poincare_maps)
    call unlink(filenames%prism_moments)
    call unlink(filenames%prism_moments_summed_squares)
    call unlink(filenames%vertex_coordinates)
    call unlink(filenames%vertex_indices)
    call unlink(filenames%prism_volumes)
    call unlink(filenames%fourier_moments)
    call unlink(filenames%refined_prism_volumes)
    call unlink(filenames%electric_potential)
    call unlink(filenames%boltzmann_density)
    call unlink(filenames%divertor_intersections)
    call unlink(filenames%tetr_moments)

end subroutine unlink_files

end module boltzmann_writing_data_mod