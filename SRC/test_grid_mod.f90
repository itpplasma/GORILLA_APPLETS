!
module test_grid_mod
!
    double precision, dimension(:,:),allocatable :: knots_EIRENE
    integer, dimension(:,:), allocatable :: triangles_EIRENE,triangles_EIRENE_region,triangle_back_interp_EIRENE, &
                                            & quadrangles_EIRENE_region
!
    logical,dimension(:),allocatable :: boole_region
    integer :: num_triangles_region
!
    contains
!
    subroutine test_grid()
!
        use orbit_timestep_gorilla_mod, only: initialize_gorilla,check_coordinate_domain
        use tetra_grid_settings_mod, only: n_field_periods
        use constants, only: ev2erg,pi
        use tetra_physics_mod, only: particle_mass
        use pusher_tetra_rk_mod, only: find_tetra
!
        implicit none
!
        integer :: i,j,k
        double precision :: s_0,theta_0,phi_0,energy_eV_0,pitchpar_0
        double precision, dimension(3) :: x
        double precision :: vmod, vpar, vperp
        integer :: ind_tetr,iface
!
        integer :: n_points
        double precision,dimension(:),allocatable :: bmod_vec,alpha_0_vec
!
        !Initialize GORILLA
        call initialize_gorilla()
!
        print *, 'Number of field periods', n_field_periods
!
        !Compute velocity module from kinetic energy dependent on particle species
        energy_eV_0 = 3.d3
        vmod=sqrt(2.d0*energy_eV_0*ev2erg/particle_mass)
!
        s_0 = 0.5d0
        theta_0 = 0.d0
        phi_0 = 0.d0
!
        !--- Find tetrahedron for starting positions by neglecting electrostatic potential energy
!
        pitchpar_0 = 0.8d0
        vpar = pitchpar_0*vmod
        vperp = sqrt((1.d0-pitchpar_0**2)*vmod**2)
!
        !Scan grid
        n_points = 1000
        allocate(bmod_vec(n_points+1),alpha_0_vec(n_points+1))
!
        alpha_0_vec(:) = [(2.d0*pi/dble(n_points)*dble(i), i = 0,n_points)]
!
        do i = 1,n_points+1
!
            !Define start position
            x(1) = 0.5d0+0.25d0*cos(alpha_0_vec(i))
            x(2) = sin(alpha_0_vec(i))
            x(3) = 2.d0*pi !phi_0
!
            !Check coordinate domain (optionally perform modulo operation)
            call check_coordinate_domain(x)
!
            !Find tetrahedron index and face index for position x
            call find_tetra(x,vpar,vperp,ind_tetr,iface)
!
            bmod_vec(i) = bmod_func(x,ind_tetr)
!
        enddo
!
        do i = 1,n_points+1
            write(100,*) bmod_vec(i)
        enddo
!
    end subroutine test_grid
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function bmod_func(x,ind_tetr)
!
        use tetra_physics_mod, only: tetra_physics
!
        implicit none
!
        double precision :: bmod_func
        integer, intent(in) :: ind_tetr
        double precision, dimension(3),intent(in) :: x
        double precision, dimension(3) :: z
!
        z = x-tetra_physics(ind_tetr)%x1
        bmod_func = tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*z)
!
    end function bmod_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine test_west()
!
        use orbit_timestep_gorilla_mod, only: initialize_gorilla
        use tetra_grid_settings_mod, only: n_field_periods
        use tetra_physics_mod, only: tetra_physics
        use supporting_functions_mod, only: sym_flux_in_cyl
        use tetra_grid_mod, only: verts_rphiz
!
        implicit none
!
        double precision :: max_dt_dtau_const,min_dt_dtau_const,max_bmod,min_bmod,min_R,max_R,min_Z,max_Z
        integer :: i,n_start
!
        !Initialize GORILLA
        call initialize_gorilla()
!
        print *, 'Number of field periods', n_field_periods
!
        !Min and Max of dt_dtau_const
        min_dt_dtau_const = minval(tetra_physics(:)%dt_dtau_const)
        max_dt_dtau_const = maxval(tetra_physics(:)%dt_dtau_const)
!
        print *, 'min_dt_dtau_const',min_dt_dtau_const
        print *, 'max_dt_dtau_const',max_dt_dtau_const
!
        !Min and Max of bmod
        min_bmod = minval(tetra_physics(:)%bmod1)
        max_bmod = maxval(tetra_physics(:)%bmod1)
!
        print *, 'min_bmod',min_bmod
        print *, 'max_bmod',max_bmod
!
!        call sym_flux_in_cyl('orbit_start_sthetaphilambda.dat', &
!            & 'orbit_start_rphizlambda.dat',1)

        !Min and Max of R and Z
        min_R = minval(verts_rphiz(1,:))
        max_R = maxval(verts_rphiz(1,:))
        min_Z = minval(verts_rphiz(3,:))
        max_Z = maxval(verts_rphiz(3,:))
!
        print *, 'min_R',min_R
        print *, 'max_R',max_R
        print *, 'min_Z',min_Z
        print *, 'max_Z',max_Z
!
        !Write file with starting positions
        n_start = 100
!
        open(unit=299, file='orbit_start_rphizlambda.dat', status='unknown')
        do i = 1,n_start
            write(299,*) min_R + 10.d0 + (max_R-min_R-20.d0)/dble(n_start)*dble(i), 0.1d0, 0.1d0, 0.7d0
        enddo
        close(299)
!
    end subroutine test_west
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine test_eirene_grid
!
        use orbit_timestep_gorilla_mod, only: initialize_gorilla
        use supporting_functions_mod, only: logical2integer
!
        use tetra_grid_mod, only: Rmin,Rmax,Zmin,Zmax
!
        implicit none
!
        integer :: i,j,k,i_region,t_type
        integer :: file_id_knots,file_id_triangles,file_id_triangle_back_interp
        integer, dimension(2) :: shape_knots,shape_triangles,shape_triangle_back_interp
        character(50) :: filename_knots,filename_triangles,filename_triangle_back_interp
!
        double precision :: r,z,A_r,A_phi,A_z,bmod
!
        !Initialize GORILLA
        call initialize_gorilla()
!
        !---------------------------------------------------------------------------------------------------------------------!
        !Load EIRENE mesh data
!
        !Define file ids and filenames for EIRENE mesh data
        file_id_knots = 501
        filename_knots = './EIRENE_MESH/knots.dat'
!
        file_id_triangles = 502
        filename_triangles = './EIRENE_MESH/triangles.dat'
!
        file_id_triangle_back_interp = 503
        filename_triangle_back_interp = './EIRENE_MESH/triangle_back_interp.dat'
!
        !Load knots
        open(unit=file_id_knots, file=filename_knots, status='unknown')
        read(file_id_knots,*) shape_knots
        allocate(knots_EIRENE(shape_knots(1),shape_knots(2)))
        do i = 1,shape_knots(1)
            read(file_id_knots,*) knots_EIRENE(i,:)
        enddo
        close(file_id_knots)
!print *, 'knots_EIRENE(shape_knots(1),:)',knots_EIRENE(shape_knots(1),:)
!
        !Load triangles
        open(unit=file_id_triangles, file=filename_triangles, status='unknown')
        read(file_id_triangles,*) shape_triangles
        allocate(triangles_EIRENE(shape_triangles(1),shape_triangles(2)))
        do i = 1,shape_triangles(1)
            read(file_id_triangles,*) triangles_EIRENE(i,:)
        enddo
        close(file_id_triangles)
print *, 'triangles_EIRENE(1,:)',triangles_EIRENE(1,:)
!
        !Load triangle_back_interp
        open(unit=file_id_triangle_back_interp, file=filename_triangle_back_interp, status='unknown')
        read(file_id_triangle_back_interp,*) shape_triangle_back_interp
        allocate(triangle_back_interp_EIRENE(shape_triangle_back_interp(1),shape_triangle_back_interp(2)))
        do i = 1,shape_triangle_back_interp(1)
            read(file_id_triangle_back_interp,*) triangle_back_interp_EIRENE(i,:)
        enddo
        close(file_id_triangle_back_interp)
print *, 'triangle_back_interp_EIRENE(shape_triangle_back_interp(1))',triangle_back_interp_EIRENE(shape_triangle_back_interp(1),1)
!
        !---------------------------------------------------------------------------------------------------------------------!
        !Select region
!
        allocate(boole_region(shape_triangle_back_interp(1)))
!

!Loop over regions
do j = 1,maxval(triangle_back_interp_EIRENE(:,1))
        i_region = j
!
print *, 'Region: ',i_region
!
        where (triangle_back_interp_EIRENE(:,1) == i_region)
            boole_region = .true.
        else where
            boole_region = .false.
        end where
!
        num_triangles_region = sum(logical2integer(boole_region))
        allocate(triangles_EIRENE_region(num_triangles_region,shape_triangles(2)))
        allocate(quadrangles_EIRENE_region(num_triangles_region,2))
        k = 1
        do i = 1,shape_triangles(1)
            if(boole_region(i)) then
                triangles_EIRENE_region(k,:) = triangles_EIRENE(i,:)
                quadrangles_EIRENE_region(k,:) = triangle_back_interp_EIRENE(i,[2,3])
                k = k+1
            endif
        enddo
!
print *, 'Number of triangles in region: ',num_triangles_region
!
        call test_triangle_type_region()
!
        deallocate(triangles_EIRENE_region,quadrangles_EIRENE_region)

enddo !Loop over regions

        !---------------------------------------------------------------------------------------------------------------------!
        !Obtain poloidal flux
!
!do j = 1,10 !num_triangles_region
!print *, 'j = ',j,'quadrangles_EIRENE_region(j,:)',quadrangles_EIRENE_region(j,:)
!        do i = 1,3
!            r = knots_EIRENE(triangles_EIRENE_region(j,i),1)
!            z = knots_EIRENE(triangles_EIRENE_region(j,i),2)
!!
!            call vector_potential_rz(r,z,A_phi)
!!
!print *, 'i',i,'A_phi',A_phi,'r',r,'z',z
!!write(777,*) r,z,A_phi
!        enddo
!print *, 'now'
!t_type =triangle_type(j)
!print *,'triangle_type',t_type
!print *, ''
!enddo

!call test_triangle_type_region()

!Write out ASDEX data
!do i = 1,100
!    R = Rmin + (Rmax-Rmin)/dble(100)*dble(i)
!    do j = 1,100
!        Z = Zmin + (Zmax-Zmin)/dble(100)*dble(j)
!        call vector_potential_rz(r,z,A_r,A_phi,A_z,bmod)
!        write(777,*) r,z,A_phi
!    enddo
!enddo
!
    end subroutine test_eirene_grid
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine vector_potential_rz(r,z,A_phi)
!
        use field_eq_mod, only : rtf,btf,psif
!
        implicit none
!
        double precision,intent(in) :: r,z
        double precision, intent(out) :: A_phi
!        double precision, intent(out) :: A_r,A_z,bmod
        double precision :: phi
        double precision :: B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ    &
                            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
!
        phi = 0.d0
!
        call field(r,phi,z,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ  &
                  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
        !bmod=sqrt(B_r**2+B_p**2+B_z**2)
!
        !A_r=0.d0
        A_phi=psif
        !A_z=-rtf*btf*log(r)
!
    end subroutine vector_potential_rz
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function triangle_type(i_triangle_EIRENE_region)
!
        implicit none
!
        integer :: triangle_type,j,i_min_diff
        integer, intent(in) :: i_triangle_EIRENE_region
        double precision :: r,z
        double precision,dimension(3) :: A_phi_vec,diff_Aphi_vec
!
        do j = 1,3
            r = knots_EIRENE(triangles_EIRENE_region(i_triangle_EIRENE_region,j),1)
            z = knots_EIRENE(triangles_EIRENE_region(i_triangle_EIRENE_region,j),2)
    !
            call vector_potential_rz(r,z,A_phi_vec(j))
        enddo
!
        diff_Aphi_vec(1) = abs(A_phi_vec(1)-A_phi_vec(2))
        diff_Aphi_vec(2) = abs(A_phi_vec(2)-A_phi_vec(3))
        diff_Aphi_vec(3) = abs(A_phi_vec(3)-A_phi_vec(1))
!
        i_min_diff = minloc(diff_Aphi_vec,1)
!
        select case(i_min_diff)
            case(1)
                if( abs(A_phi_vec(1)).gt.abs(A_phi_vec(3))) then
                    triangle_type = 1
                else
                    triangle_type = 2
                endif
            case(2)
                if( abs(A_phi_vec(2)).gt.abs(A_phi_vec(1))) then
                    triangle_type = 1
                else
                    triangle_type = 2
                endif
            case(3)
                if( abs(A_phi_vec(3)).gt.abs(A_phi_vec(2))) then
                    triangle_type = 1
                else
                    triangle_type = 2
                endif
        end select
!
    end function triangle_type
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine test_triangle_type_region()
!
        implicit none
!
        integer :: i,j,k,counter_quadrangles,counter_different_triangle_types
        double precision, dimension(2) :: cur_quadrangle
        double precision :: r,z,A_phi
!
        counter_quadrangles = 0
        counter_different_triangle_types = 0
!
        do i = 1,num_triangles_region-1
            cur_quadrangle = quadrangles_EIRENE_region(i,:)
            do j = i+1,num_triangles_region
                !if (i == j) cycle
                if( all(cur_quadrangle(:) == quadrangles_EIRENE_region(j,:) )) then
                    counter_quadrangles = counter_quadrangles + 1
!
!                    !Write coordinates of quadranagles
!                    write(778,*) knots_EIRENE(triangles_EIRENE_region(i,:),1), knots_EIRENE(triangles_EIRENE_region(j,:),1)
!                    write(779,*) knots_EIRENE(triangles_EIRENE_region(i,:),2), knots_EIRENE(triangles_EIRENE_region(j,:),2)
!
                    !Check, if triangle types differ within one quadrangle
                    if ( triangle_type(i).ne.triangle_type(j) ) then
                        counter_different_triangle_types = counter_different_triangle_types + 1
                    else
!
!Write coordinates of quadranagles
write(750,*) knots_EIRENE(triangles_EIRENE_region(i,:),1)
write(750,*) knots_EIRENE(triangles_EIRENE_region(j,:),1)
write(751,*) knots_EIRENE(triangles_EIRENE_region(i,:),2)
write(751,*) knots_EIRENE(triangles_EIRENE_region(j,:),2)
!
!                        print *, 'ERROR: Found quadrangle, where both triangles are of the same type.'
!                        print *, 'i = ',i,'quadrangles_EIRENE_region(i,:)',quadrangles_EIRENE_region(i,:)
!                        do k = 1,3
!                            r = knots_EIRENE(triangles_EIRENE_region(i,k),1)
!                            z = knots_EIRENE(triangles_EIRENE_region(i,k),2)
!                            call vector_potential_rz(r,z,A_phi)
!                            print *, 'k',k,'A_phi',A_phi,'r',r,'z',z
!                        enddo
!                        print *,'triangle_type',triangle_type(i)
!                        print *,''
!!
!                        print *, 'j = ',j,'quadrangles_EIRENE_region(j,:)',quadrangles_EIRENE_region(j,:)
!                        do k = 1,3
!                            r = knots_EIRENE(triangles_EIRENE_region(j,k),1)
!                            z = knots_EIRENE(triangles_EIRENE_region(j,k),2)
!                            call vector_potential_rz(r,z,A_phi)
!                            print *, 'k',k,'A_phi',A_phi,'r',r,'z',z
!                        enddo
!                        print *,'triangle_type',triangle_type(j)
!                        print *,''
!                        print *,''
                    endif
!
                endif !Quadrangle
            enddo
        enddo
!
!        print *, 'number of triangles', num_triangles_region
!        print *, 'number of quadrangles', counter_quadrangles
!        print *, 'number of different triangle types', counter_different_triangle_types
!
        if(counter_quadrangles.eq.counter_different_triangle_types) then
            print *, 'Test of triangle types in selected region was successful.'
            print *, ''
        else
            print *, 'Error: Test of triangle types in selected region was NOT successful.'
            print *, 'Number of non correct triangle types', counter_quadrangles -counter_different_triangle_types, '/', &
                    & counter_quadrangles
            print *, ''
!            stop
        endif
!
    end subroutine test_triangle_type_region
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function check_orientation(i_triangle)
        implicit none
        integer, intent(in) :: i_triangle
        double precision :: a,b,c,d, angle
        logical :: check_orientation

        a = knots_EIRENE(triangles_EIRENE(i_triangle,2),1) - knots_EIRENE(triangles_EIRENE(i_triangle,1),1)
        b = knots_EIRENE(triangles_EIRENE(i_triangle,2),2) - knots_EIRENE(triangles_EIRENE(i_triangle,1),2)
        c = knots_EIRENE(triangles_EIRENE(i_triangle,3),1) - knots_EIRENE(triangles_EIRENE(i_triangle,1),1)
        d = knots_EIRENE(triangles_EIRENE(i_triangle,3),2) - knots_EIRENE(triangles_EIRENE(i_triangle,1),2)

        ! perpendicular dot product -> sign of shortest angle between two vectors
        angle = a*d - b*c
        check_orientation = angle > 0 ! -> then counterclockwise (mathematical positiv)

    end function check_orientation
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module test_grid_mod
