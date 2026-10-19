!Radial electrostatic potential Phi0(s) on a flux-coordinate mesh.
!
!On flux-coordinate meshes (coord_system 2, x1 = s; e.g. the Boozer chartmap,
!grid_kind 6) GORILLA's default potential is phi = eps_Phi*A_theta, i.e. a
!constant dPhi/ds. This module replaces the vertex potential by a prescribed
!radial profile Phi0(s) [statV] and rebuilds the mesh physics, so the ExB
!drift of the orbits follows the toroidal rotation frequency
!    Omega_E(s) = c*dPhi0/dA_phi
!(field-line label alpha = theta - iota*phi changes at -iota*Omega_E).
!Within a tetrahedron the potential is linear, so dPhi0/ds is the secant
!slope of the radial cell. The caller passes the same Phi0(s) that enters
!its background distribution, e.g. the profile potential of profile_data_mod.
module flux_mesh_potential_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    abstract interface
        function radial_potential(s) result(phi)
            import :: dp
            real(dp), intent(in) :: s
            real(dp) :: phi
        end function radial_potential
    end interface

    public :: radial_potential, set_flux_mesh_potential

contains

    subroutine set_flux_mesh_potential(phi0, ipert, bmod_multiplier)

        use gorilla_settings_mod, only: coord_system
        use tetra_grid_mod, only: nvert, verts_sthetaphi
        use tetra_physics_mod, only: phi_elec, make_tetra_physics
        use tetra_physics_poly_precomp_mod, only: make_precomp_poly

        procedure(radial_potential) :: phi0
        integer, intent(in) :: ipert
        real(dp), intent(in), optional :: bmod_multiplier

        integer :: iv

        if (coord_system /= 2) error stop &
            'set_flux_mesh_potential: needs flux coordinates (coord_system = 2)'
        if (.not. allocated(phi_elec)) error stop &
            'set_flux_mesh_potential: initialize_gorilla must run first'

        do iv = 1, nvert
            phi_elec(iv) = phi0(verts_sthetaphi(1, iv))
        end do
        print '(a,2es14.6)', ' Mesh potential Phi0(s) on the vertices, min/max [statV]:', &
            minval(phi_elec), maxval(phi_elec)

        if (present(bmod_multiplier)) then
            call make_tetra_physics(coord_system, ipert, bmod_multiplier, &
                                    boole_keep_phi_elec=.true.)
        else
            call make_tetra_physics(coord_system, ipert, boole_keep_phi_elec=.true.)
        end if
        call make_precomp_poly()

    end subroutine set_flux_mesh_potential

end module flux_mesh_potential_mod
