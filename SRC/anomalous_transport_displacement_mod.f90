module anomalous_transport_displacement_mod
!
! Module for anomalous transport displacement calculations.
!
! Contains routines for:
!   - Straight-line radial displacement for anomalous diffusion
!   - Diffusion tensor Cholesky decomposition
!   - Correction velocity computation
!   - Tetrahedron face intersection finding
!

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine carry_out_anomalous_transport_displacement(x, ind_tetr, iface, dt)
!
! Computes the displacement vector for anomalous transport and applies it.
!
! This subroutine calculates the perpendicular diffusion tensor and then
! calls displace_by_straight_line to perform the actual displacement
! through the tetrahedral mesh.
!
! The displacement follows:
!   Delta x^i = sqrt(2 * dt) * alpha^{ij} * xi_j + V_c^i * dt
!
! where x^1 = R, x^2 = phi, x^3 = Z (cylindrical coordinates),
! alpha is the Cholesky decomposition of D_perp, xi_j are random numbers
! with zero mean and unit variance, and V_c is the correction velocity.
!
! Input/Output:
!   x(3)      - Particle position in (R, phi, Z) coordinates
!   ind_tetr  - Current tetrahedron index (updated after displacement)
!   iface     - Current face index (updated after displacement)
!   dt        - Time step
!
    use tetra_physics_mod, only: tetra_physics
    use collis_ions, only: getran
    use gorilla_applets_types_mod, only: g, counter_t
    use tetra_grid_settings_mod, only: grid_kind
    use utils_orbit_timestep_mod, only: identify_particles_entering_annulus
    use find_tetra_mod, only: find_tetra

    real(dp), dimension(3), intent(inout) :: x
    integer, intent(inout)                :: ind_tetr, iface
    real(dp), intent(in)                  :: dt

    real(dp), dimension(3) :: displacement, xi, V_c, x_new
    real(dp), parameter :: diffusion_coefficient = 1.0d-4  ! TODO: make this an input parameter

    ! Cholesky decomposition of the perpendicular diffusion tensor
    real(dp), dimension(3,3) :: alpha_perp_mat
    real(dp), dimension(3) :: h_contra, x_local
    real(dp) :: R_local, sqrt_2dt
    logical :: boole_lost_inside
    type(counter_t) :: dummy_counter
    real(dp) :: vpar_dummy, vperp_dummy

    ! Compute position relative to first vertex
    x_local = x - tetra_physics(ind_tetr)%x1

    ! Interpolate R (major radius) at current position
    R_local = tetra_physics(ind_tetr)%R1 + dot_product(tetra_physics(ind_tetr)%gR, x_local)

    ! Interpolate h (unit vector along B) at current position using linear approximation
    ! h_i(x) = h_i_1 + grad(h_i) . (x - x1)
    h_contra(1) =  tetra_physics(ind_tetr)%h1_1 + dot_product(tetra_physics(ind_tetr)%gh1, x_local)
    h_contra(2) = (tetra_physics(ind_tetr)%h2_1 + dot_product(tetra_physics(ind_tetr)%gh2, x_local)) / (R_local**2)
    h_contra(3) =  tetra_physics(ind_tetr)%h3_1 + dot_product(tetra_physics(ind_tetr)%gh3, x_local)

    ! Compute Cholesky decomposition of the diffusion tensor
    call compute_diffusion_cholesky(h_contra, R_local, diffusion_coefficient, alpha_perp_mat)

    ! Compute correction velocity
    call compute_correction_velocity(ind_tetr, h_contra, R_local, diffusion_coefficient, V_c)

    ! Generate random numbers with zero mean and unit variance
    call getran(0, xi(1))
    call getran(0, xi(2))
    call getran(0, xi(3))

    ! Compute displacement: Delta x^i = sqrt(2 * dt) * alpha^{ij} * xi_j + V_c^i * dt
    ! Note: alpha_perp_mat is lower triangular, so alpha^{ij} * xi_j expands as:
    ! displacement(1) = alpha(1,1)*xi(1)
    ! displacement(2) = alpha(2,1)*xi(1) + alpha(2,2)*xi(2)
    ! displacement(3) = alpha(3,1)*xi(1) + alpha(3,2)*xi(2) + alpha(3,3)*xi(3)
    sqrt_2dt = sqrt(2.0_dp * dt)
    displacement(1) = sqrt_2dt * alpha_perp_mat(1,1) * xi(1) + V_c(1) * dt
    displacement(2) = sqrt_2dt * (alpha_perp_mat(2,1) * xi(1) + alpha_perp_mat(2,2) * xi(2)) + V_c(2) * dt
    displacement(3) = sqrt_2dt * (alpha_perp_mat(3,1) * xi(1) + alpha_perp_mat(3,2) * xi(2) &
                                + alpha_perp_mat(3,3) * xi(3)) + V_c(3) * dt

    call displace_by_straight_line(x, ind_tetr, iface, displacement)

    ! Check if particle was lost and attempt to push across the annulus if applicable
    if (ind_tetr.eq.-1) then
        if ((grid_kind.eq.2).or.(grid_kind.eq.3)) then
            ! Initialize dummy counter (not used in this context)
            dummy_counter%lost_inside = 0
            call identify_particles_entering_annulus(x, dummy_counter, boole_lost_inside)
            if (boole_lost_inside) then
                ! Reflect particle across the magnetic axis
                x_new = 3*(/g%raxis, x(2), g%zaxis/) - 2*x
                vpar_dummy = 0.0_dp
                vperp_dummy = 0.0_dp
                call find_tetra(x_new, vpar_dummy, vperp_dummy, ind_tetr, iface)
                if (ind_tetr.ne.-1) then
                    x = x_new
                    print*, "particle pushing across the hole during anomalous transport displacement was successful"
                endif
                ! If ind_tetr is still -1, particle remains lost
            endif
        endif
    endif

end subroutine carry_out_anomalous_transport_displacement

! ====================================================================
subroutine compute_diffusion_cholesky(h_contra, R_local, D_perp, alpha_perp_mat)
!
! Computes the Cholesky decomposition of the perpendicular diffusion tensor.
!
! The perpendicular diffusion tensor is:
!   D_perp^{ij} = D_perp * (g^{ij} - h^i * h^j)
!
! where g^{ij} is the contravariant metric tensor (cylindrical coordinates:
! g^{11} = g^{33} = 1, g^{22} = 1/R^2, off-diagonal = 0) and h^i are the
! contravariant components of the unit vector in the magnetic field direction.
!
! The Cholesky decomposition produces a lower triangular matrix alpha such that
! D_perp = alpha * alpha^T.
!
! Input:
!   h_contra(3) - Contravariant components of the magnetic field unit vector
!   R_local     - Local major radius
!   D_perp      - Perpendicular diffusion coefficient (scalar)
!
! Output:
!   alpha_perp_mat(3,3) - Lower triangular Cholesky factor
!
    real(dp), dimension(3), intent(in)    :: h_contra
    real(dp), intent(in)                  :: R_local, D_perp
    real(dp), dimension(3,3), intent(out) :: alpha_perp_mat

    real(dp), dimension(3,3) :: D_perp_mat, g_mat
    real(dp) :: arg_sqrt, eps_cholesky

    ! Small regularization parameter for Cholesky decomposition
    eps_cholesky = 1.0d-20

    ! Contravariant metric tensor for cylindrical coordinates (R, phi, Z)
    ! g^{11} = 1, g^{22} = 1/R^2, g^{33} = 1, off-diagonal = 0
    g_mat = 0.0_dp
    g_mat(1,1) = 1.0_dp
    g_mat(2,2) = 1.0_dp / (R_local**2)
    g_mat(3,3) = 1.0_dp

    ! Perpendicular diffusion tensor: D_perp^{ij} = D_perp * (g^{ij} - h^i * h^j)
    D_perp_mat(1,1) = D_perp * (g_mat(1,1) - h_contra(1) * h_contra(1))
    D_perp_mat(1,2) = D_perp * (g_mat(1,2) - h_contra(1) * h_contra(2))
    D_perp_mat(1,3) = D_perp * (g_mat(1,3) - h_contra(1) * h_contra(3))
    D_perp_mat(2,1) = D_perp * (g_mat(2,1) - h_contra(2) * h_contra(1))
    D_perp_mat(2,2) = D_perp * (g_mat(2,2) - h_contra(2) * h_contra(2))
    D_perp_mat(2,3) = D_perp * (g_mat(2,3) - h_contra(2) * h_contra(3))
    D_perp_mat(3,1) = D_perp * (g_mat(3,1) - h_contra(3) * h_contra(1))
    D_perp_mat(3,2) = D_perp * (g_mat(3,2) - h_contra(3) * h_contra(2))
    D_perp_mat(3,3) = D_perp * (g_mat(3,3) - h_contra(3) * h_contra(3))

    ! Cholesky decomposition: alpha_perp such that D_perp = alpha_perp * alpha_perp^T
    ! Row 1
    arg_sqrt = max(D_perp_mat(1,1), eps_cholesky)
    alpha_perp_mat(1,1) = sqrt(arg_sqrt)
    alpha_perp_mat(1,2) = 0.0_dp
    alpha_perp_mat(1,3) = 0.0_dp
    ! Row 2
    alpha_perp_mat(2,1) = D_perp_mat(1,2) / alpha_perp_mat(1,1)
    arg_sqrt = max(D_perp_mat(2,2) - alpha_perp_mat(2,1)**2, eps_cholesky)
    alpha_perp_mat(2,2) = sqrt(arg_sqrt)
    alpha_perp_mat(2,3) = 0.0_dp
    ! Row 3
    alpha_perp_mat(3,1) = D_perp_mat(1,3) / alpha_perp_mat(1,1)
    alpha_perp_mat(3,2) = (D_perp_mat(2,3) - alpha_perp_mat(2,1) * alpha_perp_mat(3,1)) / alpha_perp_mat(2,2)
    arg_sqrt = max(D_perp_mat(3,3) - alpha_perp_mat(3,1)**2 - alpha_perp_mat(3,2)**2, eps_cholesky)
    alpha_perp_mat(3,3) = sqrt(arg_sqrt)

end subroutine compute_diffusion_cholesky

! ====================================================================
subroutine compute_correction_velocity(ind_tetr, h_contra, R_local, D_perp, V_c)
!
! Computes the correction velocity for the stochastic differential equation.
!
! The correction velocity is:
!   V_c^i = V^i + (1/sqrt(g)) * d/dx^j (sqrt(g) * D_perp^{ij})
!         = V^i + dD_perp^{ij}/dx^j + D_perp^{i1}/R
!
! where V^i is the drift velocity (set to zero for now), and the second term
! comes from the coordinate system (cylindrical: sqrt(g) = R).
!
! The derivative of the diffusion tensor is:
!   dD_perp^{ij}/dx^j = -D_perp * (h^i * dh^j/dx^j + h^j * dh^i/dx^j)
!
! Input:
!   ind_tetr    - Current tetrahedron index
!   h_contra(3) - Contravariant components of the magnetic field unit vector
!   R_local     - Local major radius
!   D_perp      - Perpendicular diffusion coefficient (scalar)
!
! Output:
!   V_c(3)      - Correction velocity
!
    use tetra_physics_mod, only: tetra_physics

    integer, intent(in)                   :: ind_tetr
    real(dp), dimension(3), intent(in)    :: h_contra
    real(dp), intent(in)                  :: R_local, D_perp
    real(dp), dimension(3), intent(out)   :: V_c

    real(dp), dimension(3) :: V, partial_D_perp_ij_partial_xj, D_perp_i1
    real(dp), dimension(3,3) :: grad_h_contra  ! grad_h_contra(i,j) = dh^i/dx^j
    real(dp) :: div_h_contra

    ! Drift velocity (set to zero for now)
    V = 0.0_dp

    ! Gradient of contravariant h components: grad_h_contra(i,j) = dh^i/dx^j
    ! For cylindrical coordinates, h^i = g^{ij} h_j, so:
    ! h^1 = h_1, h^2 = h_2/R^2, h^3 = h_3
    ! dh^1/dx^j = dh_1/dx^j = gh1(j)
    ! dh^2/dx^j = d(h_2/R^2)/dx^j = gh2(j)/R^2 - 2*h^2/R * dR/dx^j
    ! dh^3/dx^j = dh_3/dx^j = gh3(j)
    grad_h_contra(1,:) = tetra_physics(ind_tetr)%gh1
    grad_h_contra(2,:) = tetra_physics(ind_tetr)%gh2 / (R_local**2) &
                       - 2.0_dp * h_contra(2) / R_local * tetra_physics(ind_tetr)%gR
    grad_h_contra(3,:) = tetra_physics(ind_tetr)%gh3

    ! Divergence of h^j: div_h_contra = dh^j/dx^j (sum over j)
    div_h_contra = grad_h_contra(1,1) + grad_h_contra(2,2) + grad_h_contra(3,3)

    ! Derivative of perpendicular diffusion tensor: dD_perp^{ij}/dx^j = -D_perp * (h^i * dh^j/dx^j + h^j * dh^i/dx^j)
    ! Summed over j, this gives a vector (one component for each i)
    partial_D_perp_ij_partial_xj(1) = -D_perp * ( &
        h_contra(1) * div_h_contra + &
        h_contra(1) * grad_h_contra(1,1) + h_contra(2) * grad_h_contra(1,2) + h_contra(3) * grad_h_contra(1,3))
    partial_D_perp_ij_partial_xj(2) = -D_perp * ( &
        h_contra(2) * div_h_contra + &
        h_contra(1) * grad_h_contra(2,1) + h_contra(2) * grad_h_contra(2,2) + h_contra(3) * grad_h_contra(2,3))
    partial_D_perp_ij_partial_xj(3) = -D_perp * ( &
        h_contra(3) * div_h_contra + &
        h_contra(1) * grad_h_contra(3,1) + h_contra(2) * grad_h_contra(3,2) + h_contra(3) * grad_h_contra(3,3))

    ! Compute D_perp^{i1} = D_perp * (g^{i1} - h^i * h^1)
    ! For cylindrical coordinates: g^{11} = 1, g^{21} = g^{31} = 0
    ! So: D_perp^{11} = D_perp * (1 - h^1 * h^1)
    !     D_perp^{21} = D_perp * (0 - h^2 * h^1) = -D_perp * h^2 * h^1
    !     D_perp^{31} = D_perp * (0 - h^3 * h^1) = -D_perp * h^3 * h^1
    D_perp_i1(1) = D_perp * (1.0_dp - h_contra(1) * h_contra(1))
    D_perp_i1(2) = -D_perp * h_contra(2) * h_contra(1)
    D_perp_i1(3) = -D_perp * h_contra(3) * h_contra(1)

    ! Correction velocity: V_c^i = V^i + dD_perp^{ij}/dx^j + D_perp^{i1}/R
    V_c(1) = V(1) + partial_D_perp_ij_partial_xj(1) + D_perp_i1(1) / R_local
    V_c(2) = V(2) + partial_D_perp_ij_partial_xj(2) + D_perp_i1(2) / R_local
    V_c(3) = V(3) + partial_D_perp_ij_partial_xj(3) + D_perp_i1(3) / R_local

end subroutine compute_correction_velocity

! ====================================================================
subroutine displace_by_straight_line(x, ind_tetr, iface, displacement)
!
! Displaces a particle along a straight line through the tetrahedral mesh.
!
! This version traverses tetrahedra along the displacement path instead of
! using find_tetra after each displacement. This is much more efficient as
! it only needs to check faces of the current tetrahedron rather than
! searching through all tetrahedra.
!
! Input/Output:
!   x(3)        - Particle position in cylindrical (R, phi, Z) coordinates
!   ind_tetr    - Current tetrahedron index (updated after displacement)
!   iface       - Current face index (updated after displacement)
!   displacement(3) - Displacement vector in cylindrical (R, phi, Z) coordinates
!
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_mod, only: tetra_grid
    use pusher_tetra_func_mod, only: pusher_handover2neighbour
    use constants, only: eps

    real(dp), dimension(3), intent(inout) :: x
    integer, intent(inout)                :: ind_tetr, iface
    real(dp), dimension(3), intent(in)    :: displacement

    real(dp) :: remaining_distance, distance_to_face, dist_min, total_distance
    real(dp), dimension(3) :: direction, z_rel
    integer :: hit_face, iper_phi, ind_tetr_old
    integer :: max_crossings, crossing_count

    ! Compute total displacement distance and direction
    total_distance = sqrt(displacement(1)**2 + displacement(2)**2 + displacement(3)**2)
    if (total_distance < 1.0d-15) return  ! No displacement

    direction = displacement / total_distance
    remaining_distance = total_distance

    ! Maximum number of tetrahedron crossings to prevent infinite loops
    max_crossings = 10000
    crossing_count = 0

    ! Use relative tolerance based on tetrahedron size, consistent with find_tetra_face_intersection
    dist_min = eps * abs(tetra_physics(ind_tetr)%dist_ref)

    ! Traverse tetrahedra until displacement is complete
    do while (remaining_distance > dist_min .and. crossing_count < max_crossings)
        crossing_count = crossing_count + 1

        ! Compute position relative to first vertex of current tetrahedron
        z_rel = x - tetra_physics(ind_tetr)%x1

        ! Find intersection with tetrahedron faces
        call find_tetra_face_intersection(z_rel, direction, ind_tetr, distance_to_face, hit_face)

        if (hit_face == 0) then
            ! No valid intersection found - this shouldn't happen if particle is inside
            ! Move the full remaining distance and exit
            x = x + remaining_distance * direction
            exit
        endif

        if (distance_to_face >= remaining_distance) then
            ! Displacement ends within this tetrahedron
            x = x + remaining_distance * direction
            remaining_distance = 0.0_dp
        else
            ! Move to the face and continue to next tetrahedron
            x = x + distance_to_face * direction

            ! Get the neighboring tetrahedron
            ind_tetr_old = ind_tetr
            call pusher_handover2neighbour(ind_tetr_old, ind_tetr, hit_face, x, iper_phi)

            ! Update iface to the entry face of the new tetrahedron
            iface = hit_face

            ! Check if we've left the domain
            if (ind_tetr == -1) then
                ! Particle has left the tetrahedral mesh
                exit
            endif

            remaining_distance = remaining_distance - distance_to_face
        endif
    enddo

end subroutine displace_by_straight_line

! ====================================================================
subroutine find_tetra_face_intersection(start_point, direction, ind_tetr, &
                                        distance_to_face, hit_face)
!
! Finds the intersection of a ray with a tetrahedron's faces.
!
! Uses the pre-computed face normals from tetra_physics to find which face
! the ray exits through. This is analogous to how the particle pusher works.
!
! The tetrahedron face normals point INWARD. A particle inside the tetrahedron
! will cross face i when the normal distance to face i becomes negative.
! For a ray from start_point in direction 'direction', we find when:
!   normal_distance(t) = normal_distance_0 + t * (anorm . direction) = 0
!
! Special handling for particles on a face:
!   When a particle is exactly on a face (normal_dist_0 ≈ 0), we need to
!   determine if it should exit through that face. We accept t_intersect >= 0
!   but when multiple faces have t ≈ 0, we pick the one we're moving towards
!   most strongly (largest negative normal_dot_dir).
!
! Face indexing in GORILLA:
!   Face 1: opposite to vertex 1 (contains vertices 2,3,4)
!   Face 2: opposite to vertex 2 (contains vertices 1,3,4)
!   Face 3: opposite to vertex 3 (contains vertices 1,2,4)
!   Face 4: opposite to vertex 4 (contains vertices 1,2,3)
!
! Input:
!   start_point(3)    - Starting position inside the tetrahedron (relative to x1)
!   direction(3)      - Direction vector of the ray (NOT necessarily unit)
!   ind_tetr          - Index of the current tetrahedron
!
! Output:
!   distance_to_face  - Distance along direction to the exit face
!                       (in units of |direction|, so actual distance = distance_to_face * |direction|)
!   hit_face          - Index of the exit face (1-4), or 0 if no valid intersection
!
    use tetra_physics_mod, only: tetra_physics
    use constants, only: eps

    real(dp), intent(in)  :: start_point(3)
    real(dp), intent(in)  :: direction(3)
    integer,  intent(in)  :: ind_tetr
    real(dp), intent(out) :: distance_to_face
    integer,  intent(out) :: hit_face

    real(dp) :: normal_dot_dir, normal_dist_0, t_intersect, dist_min
    real(dp) :: best_normal_dot_dir
    integer  :: iface

    hit_face = 0
    distance_to_face = huge(distance_to_face)
    best_normal_dot_dir = 0.0_dp

    ! Use relative tolerance based on tetrahedron size, consistent with isinside()
    dist_min = eps * abs(tetra_physics(ind_tetr)%dist_ref)

    ! Loop over all 4 faces
    do iface = 1, 4
        ! Compute: anorm . direction
        ! If negative, the ray is moving towards this face (since normals point inward)
        normal_dot_dir = dot_product(tetra_physics(ind_tetr)%anorm(:,iface), direction)

        ! Only consider faces we're moving towards (normal_dot_dir < 0 means moving towards the face)
        if (normal_dot_dir >= -dist_min) cycle

        ! Compute current normal distance to this face
        ! For face 1: includes dist_ref offset
        ! For faces 2-4: just anorm . x_local where x_local = position relative to x1
        if (iface == 1) then
            normal_dist_0 = dot_product(tetra_physics(ind_tetr)%anorm(:,1), start_point) &
                          + tetra_physics(ind_tetr)%dist_ref
        else
            normal_dist_0 = dot_product(tetra_physics(ind_tetr)%anorm(:,iface), start_point)
        endif

        ! Solve for t: normal_dist_0 + t * normal_dot_dir = 0
        ! => t = -normal_dist_0 / normal_dot_dir
        t_intersect = -normal_dist_0 / normal_dot_dir

        ! Accept intersections in forward direction (t >= -dist_min to handle numerical noise)
        if (t_intersect < -dist_min) cycle

        ! Determine if this is the best candidate
        if (t_intersect < distance_to_face - dist_min) then
            ! Clearly closer than current best
            distance_to_face = t_intersect
            hit_face = iface
            best_normal_dot_dir = normal_dot_dir
        else if (abs(t_intersect - distance_to_face) <= dist_min) then
            ! Similar distance (particle near an edge/corner) - pick face we're moving towards most strongly
            if (normal_dot_dir < best_normal_dot_dir) then
                distance_to_face = t_intersect
                hit_face = iface
                best_normal_dot_dir = normal_dot_dir
            endif
        endif
    enddo

end subroutine find_tetra_face_intersection

end module anomalous_transport_displacement_mod
