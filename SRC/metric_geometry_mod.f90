module metric_geometry_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    public :: compute_flux_coordinate_prism_volume
    public :: metric_determinant_from_local

contains

function metric_determinant_from_local(ind_tetr, z_local) result(metric_factor)

    use tetra_physics_mod, only: metric_determinant, tetra_physics

    integer, intent(in) :: ind_tetr
    real(dp), intent(in) :: z_local(3)
    real(dp) :: metric_factor
    real(dp) :: x_global(3)

    x_global = tetra_physics(ind_tetr)%x1 + z_local
    metric_factor = abs(metric_determinant(ind_tetr, x_global))

end function metric_determinant_from_local

function compute_flux_coordinate_prism_volume(ind_prism) result(prism_volume)

    use constants, only: pi
    use tetra_grid_mod, only: tetra_grid, verts_sthetaphi
    use tetra_grid_settings_mod, only: n_field_periods
    use tetra_physics_mod, only: metric_determinant, tetra_physics

    integer, intent(in) :: ind_prism
    integer :: ind_tetr
    integer :: even
    integer :: odd
    real(dp) :: basic_volume
    real(dp) :: prism_volume
    real(dp) :: verts_phi(2)
    real(dp) :: verts_s(2)
    real(dp) :: verts_theta(2)
    real(dp) :: vertex2(3)
    real(dp) :: vertex3(3)
    real(dp) :: vertex4(3)
    real(dp) :: vertex_metric_sum

    ind_tetr = (ind_prism - 1) * 3 + 1

    verts_s(1) = verts_sthetaphi(1, tetra_grid(ind_tetr)%ind_knot(1))
    verts_s(2) = verts_sthetaphi(1, tetra_grid(ind_tetr)%ind_knot(4))
    verts_theta(1) = verts_sthetaphi(2, tetra_grid(ind_tetr)%ind_knot(1))
    verts_theta(2) = verts_sthetaphi(2, tetra_grid(ind_tetr)%ind_knot(4))
    verts_phi(1) = verts_sthetaphi(3, tetra_grid(ind_tetr)%ind_knot(1))
    verts_phi(2) = verts_sthetaphi(3, tetra_grid(ind_tetr)%ind_knot(4))

    if ((verts_theta(2) - verts_theta(1)) < 0.0_dp) verts_theta(2) = verts_theta(2) + 2.0_dp * pi
    if ((verts_phi(2) - verts_phi(1)) < 0.0_dp) verts_phi(2) = verts_phi(2) + 2.0_dp * pi / real(n_field_periods, dp)

    basic_volume = (verts_s(2) - verts_s(1)) * (verts_theta(2) - verts_theta(1)) * (verts_phi(2) - verts_phi(1)) / 2.0_dp

    odd = mod(ind_prism, 2)
    even = abs(odd - 1)

    vertex2 = [ (verts_s(2) - verts_s(1)) * real(odd, dp), (verts_theta(2) - verts_theta(1)) * real(even, dp), 0.0_dp ]
    vertex3 = [ verts_s(2) - verts_s(1), verts_theta(2) - verts_theta(1), 0.0_dp ]

    vertex_metric_sum = 0.0_dp
    vertex_metric_sum = vertex_metric_sum + metric_determinant(ind_tetr, tetra_physics(ind_tetr)%x1)
    vertex_metric_sum = vertex_metric_sum + metric_determinant(ind_tetr, tetra_physics(ind_tetr)%x1 + vertex2)
    vertex_metric_sum = vertex_metric_sum + metric_determinant(ind_tetr, tetra_physics(ind_tetr)%x1 + vertex3)

    vertex2 = [ 0.0_dp, 0.0_dp, verts_phi(2) - verts_phi(1) ]
    vertex3 = [ (verts_s(2) - verts_s(1)) * real(odd, dp), (verts_theta(2) - verts_theta(1)) * real(even, dp), &
        verts_phi(2) - verts_phi(1) ]
    vertex4 = [ verts_s(2) - verts_s(1), verts_theta(2) - verts_theta(1), verts_phi(2) - verts_phi(1) ]

    vertex_metric_sum = vertex_metric_sum + metric_determinant(ind_tetr + 2, tetra_physics(ind_tetr + 2)%x1 + vertex2)
    vertex_metric_sum = vertex_metric_sum + metric_determinant(ind_tetr + 2, tetra_physics(ind_tetr + 2)%x1 + vertex3)
    vertex_metric_sum = vertex_metric_sum + metric_determinant(ind_tetr + 2, tetra_physics(ind_tetr + 2)%x1 + vertex4)

    prism_volume = basic_volume * abs(vertex_metric_sum) / 6.0_dp

end function compute_flux_coordinate_prism_volume

end module metric_geometry_mod
