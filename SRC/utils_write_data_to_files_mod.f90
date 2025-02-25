module utils_write_data_to_files_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    public :: write_data_to_files, give_file_names, unlink_files

contains

subroutine write_data_to_files

    use boltzmann_types_mod, only: filenames, in, output, moment_specs

    if (in%boole_write_vertex_indices) call write_vertex_indices

    if (in%boole_write_vertex_coordinates) call write_vertex_coordinates

    if (in%boole_write_prism_volumes) call write_prism_volumes

    if (in%boole_write_refined_prism_volumes) call write_refined_prism_volumes

    if (in%boole_write_boltzmann_density) call write_boltzmann_densities

    if (in%boole_write_electric_potential) call write_electric_potential

    if (in%boole_write_moments) then
        if (moment_specs%n_moments.gt.0) then
            call write_moments
        else
            print*, "Error: moments are not written to file because no moment was computed. Turn computation of moments on in &
                     gorilla.inp."
        endif
    endif

    if (in%boole_write_fourier_moments) then
        if (moment_specs%n_moments.gt.0) then
            call write_fourier_moments
        else
            print*, "Error: Fourier moments are not written to file because no moment was computed (and thus also no fourier &
                     moment). Turn computation of moments on in gorilla.inp."
        endif

    endif

    if (in%boole_write_exit_data) call write_exit_data

    if (in%boole_write_grid_data) call write_grid_data

end subroutine write_data_to_files

subroutine write_vertex_indices

    use tetra_grid_mod, only: ntetr, tetra_grid
    use boltzmann_types_mod, only: filenames

    integer :: vi_unit
    integer :: i

    open(newunit = vi_unit, file = filenames%vertex_indices)
    do i=1, ntetr
        write(vi_unit, *) tetra_grid(i)%ind_knot([1, 2, 3, 4])
    end do
    close(vi_unit)

end subroutine write_vertex_indices

subroutine write_vertex_coordinates

    use tetra_physics_mod, only: coord_system
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, nvert
    use boltzmann_types_mod, only: filenames

    integer :: vc_unit
    integer :: i

    open(newunit = vc_unit, file = filenames%vertex_coordinates)
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

subroutine write_prism_volumes

    use boltzmann_types_mod, only: filenames, output

    integer :: pv_unit

    open(newunit = pv_unit, file = filenames%prism_volumes)
    write(pv_unit,'(ES20.10E4)') output%prism_volumes
    close(pv_unit)

end subroutine write_prism_volumes

subroutine write_refined_prism_volumes

    use boltzmann_types_mod, only: filenames, output

    integer :: rpv_unit

    open(newunit = rpv_unit, file = filenames%refined_prism_volumes)
    write(rpv_unit,'(ES20.10E4)') output%refined_prism_volumes
    close(rpv_unit)

end subroutine write_refined_prism_volumes

subroutine write_boltzmann_densities

    use boltzmann_types_mod, only: filenames, output

    integer :: bd_unit

    open(newunit = bd_unit, file = filenames%boltzmann_density)
    write(bd_unit,'(ES20.10E4)') output%boltzmann_density
    close(bd_unit)

end subroutine write_boltzmann_densities

subroutine write_electric_potential

    use boltzmann_types_mod, only: filenames, output

    integer :: epv_unit

    open(newunit = epv_unit, file = filenames%electric_potential)
    write(epv_unit,'(ES20.10E4)') output%electric_potential
    close(epv_unit)

end subroutine write_electric_potential

subroutine write_moments

    use tetra_grid_mod, only: ntetr
    use boltzmann_types_mod, only: moment_specs, filenames, output

    integer :: p_moments_unit, pmss_unit, t_moments_unit
    integer :: l, i
    integer :: n_prisms

    open(newunit = p_moments_unit, file = filenames%prism_moments)
    open(newunit = pmss_unit, file = filenames%prism_moments_summed_squares)
    open(newunit = t_moments_unit, file = filenames%tetr_moments)

    n_prisms = ntetr/3

    if (moment_specs%n_moments.gt.0) then
        do l = 1,n_prisms
            do i = 1,moment_specs%n_moments - 1
                write(p_moments_unit,'(2ES20.10E4)',advance="no") real(output%prism_moments(i,l)), aimag(output%prism_moments(i,l))
            enddo
                write(p_moments_unit,'(2ES20.10E4)') real(output%prism_moments(moment_specs%n_moments,l)), &
                                                     aimag(output%prism_moments(moment_specs%n_moments,l))
        enddo
        if (moment_specs%boole_squared_moments) then
            do l = 1,n_prisms
                do i = 1,moment_specs%n_moments - 1
                    write(pmss_unit,'(2ES20.10E4)',advance="no") real(output%prism_moments_squared(i,l)), &
                                                                    & aimag(output%prism_moments_squared(i,l))
                enddo
                    write(pmss_unit,'(2ES20.10E4)') real(output%prism_moments_squared(moment_specs%n_moments,l)), &
                                                    aimag(output%prism_moments_squared(moment_specs%n_moments,l))
            enddo
        endif
        do l = 1,ntetr
            do i = 1,moment_specs%n_moments - 1
                write(t_moments_unit,'(2ES20.10E4)',advance="no") real(output%tetr_moments(i,l)), aimag(output%tetr_moments(i,l))
            enddo
                write(t_moments_unit,'(2ES20.10E4)') real(output%tetr_moments(moment_specs%n_moments,l)), &
                                                     aimag(output%tetr_moments(moment_specs%n_moments,l))
        enddo
    endif

    close(p_moments_unit)
    close(pmss_unit)
    close(t_moments_unit)

end subroutine write_moments

subroutine write_fourier_moments

    use tetra_grid_settings_mod, only: grid_size
    use boltzmann_types_mod, only: moment_specs, output, filenames

    integer :: fm_unit, l, i

    open(newunit = fm_unit, file = filenames%fourier_moments)
    do l = 1,moment_specs%n_triangles
        do i = 1,moment_specs%n_fourier_modes-1
            write(fm_unit,'(2ES20.10E4)',advance="no") real(output%moments_in_frequency_space(1,l,i)), &
                                                       aimag(output%moments_in_frequency_space(1,l,i))
        enddo
            write(fm_unit,'(2ES20.10E4)') real(output%moments_in_frequency_space(1,l,moment_specs%n_fourier_modes)), &
                                          aimag(output%moments_in_frequency_space(1,l,moment_specs%n_fourier_modes))
    enddo
    close(fm_unit)

end subroutine write_fourier_moments

subroutine write_exit_data

    use boltzmann_types_mod, only: exit_data, filenames, in

    integer :: ed_unit, i

    open(newunit = ed_unit, file = filenames%exit_data)
    do i = 1,in%num_particles
        write(ed_unit,*) i, exit_data%lost(i), exit_data%t_confined(i), exit_data%x(:,i), exit_data%vpar(i), &
                        exit_data%vperp(i), exit_data%integration_step(i), exit_data%phi_0_mappings(i)
    enddo
    close (ed_unit)

end subroutine write_exit_data

subroutine write_grid_data

    use tetra_grid_settings_mod, only: grid_size, n_extra_rings
    use boltzmann_types_mod, only: filenames, g
    use netcdf

    integer :: status, ncid
    integer :: dimid_2, dimid_3, varid
    
    status = nf90_create(filenames%grid_data, NF90_NETCDF4, ncid)
    call nc_check(status, 'open')

    ! define dimensions
    status = nf90_def_dim(ncid, '2', 2, dimid_2)
    status = nf90_def_dim(ncid, '3', 3, dimid_3)

    !write grid_size
    status = nf90_def_var(ncid, 'grid_size', NF90_DOUBLE, [dimid_3], varid)
    status = nf90_put_var(ncid, varid, grid_size)
    status = nf90_put_att(ncid, varid, 'comment', 'number of tetrahedra in r-, phi- and z-direction')

    !write o-point
    status = nf90_def_var(ncid, 'opoint', NF90_DOUBLE, [dimid_2], varid)
    status = nf90_put_var(ncid, varid, (/g%raxis, g%zaxis/))
    status = nf90_put_att(ncid, varid, 'comment', 'r- and z-value of o-point')

    !write minimum radius of circle centered at magnetic axis completely contained in the grid
    status = nf90_def_var(ncid, 'distance from opoint', NF90_DOUBLE, varid)
    status = nf90_put_var(ncid, varid, g%dist_from_o_point_within_grid)
    status = nf90_put_att(ncid, varid, 'comment', 'minimum radius of circle centered at &
    magnetic axis completely contained in the grid') 

    !write number of extra rings
    status = nf90_def_var(ncid, 'number of extra rings', NF90_DOUBLE, varid)
    status = nf90_put_var(ncid, varid, n_extra_rings)
    status = nf90_put_att(ncid, varid, 'comment', 'number of extra rings added close to the o-point') 

    !write range in r-direction
    status = nf90_def_var(ncid, 'r range', NF90_DOUBLE, [dimid_2], varid)
    status = nf90_put_var(ncid, varid, (/g%amin, g%amax/))
    status = nf90_put_att(ncid, varid, 'comment', 'minimum and maximum r-values') 

    !write range in z-direction
    status = nf90_def_var(ncid, 'z range', NF90_DOUBLE, [dimid_2], varid)
    status = nf90_put_var(ncid, varid, (/g%cmin, g%cmax/))
    status = nf90_put_att(ncid, varid, 'comment', 'minimum and maximum z-values') 

    status = nf90_close(ncid)
    call nc_check(status, 'close')

end subroutine write_grid_data

subroutine nc_check(status, operation)
    use netcdf, only: NF90_NOERR, nf90_strerror
    integer, intent(in) :: status
    character(len=*), intent(in) :: operation

    if (status == NF90_NOERR) return
    write (*, '("Error encountered during ", a, ": ", a)') operation, nf90_strerror(status)
    error stop
end subroutine nc_check

subroutine give_file_names

    use boltzmann_types_mod, only: filenames

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
    filenames%exit_data = 'exit_data.dat'
    filenames%grid_data = 'grid_data.nc'

end subroutine give_file_names

subroutine unlink_files

    use boltzmann_types_mod, only: filenames

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
    call unlink(filenames%exit_data)

end subroutine unlink_files

end module utils_write_data_to_files_mod