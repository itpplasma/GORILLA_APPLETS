program test_neutral_flight

! Straight-line neutral flight (neutral_flight_mod) on GORILLA's analytic circular tokamak
! (grid_kind = 5, coord_system = 1, R0 = 200 cm, a = 50 cm, no equilibrium files needed).
! Arguments: n2 (toroidal slices per field period), n_field_periods, sag tolerance for 'auto',
! expected 'auto' choice (cartesian or hybrid).
!   1. sag of the planar Cartesian faces: at least the chord sag R_max (1 - cos(dphi/2)) of the outermost
!      toroidal edge; for thin slices the twist of faces extending in R and phi dominates (~ dR dphi/4,
!      dR = R extent of a tetrahedron), so at most chord sag + dR_max dphi/2
!   2. face_hybrid wall exit vs an independent computation: the hybrid domain is the polygon of the
!      boundary vertices in (R, Z) times phi, its first exit along the line is found by bisection.
!      Launches include a toroidal line over many slices and a crossing of the field-period seam.
!   3. 'auto' choice; with face_cartesian: wall exit within a few chord sags of the hybrid one
!   4. ionisation with a constant rate nu: wall fraction exp(-nu t_wall), ionised before t_wall/2
!      1 - exp(-nu t_wall/2) (4.5 binomial sigma); the ionisation point lies in the returned
!      GORILLA tetrahedron.
!   5. recombination and ionisation (marker_species_mod): ion_to_neutral puts the particle at
!      X + (b x v_perp)/Omega with |v| kept and the displacement perpendicular to b; neutral_to_ion
!      at that point with that velocity gives back the guiding centre X and v_par, v_perp (up to the
!      variation of b over a Larmor radius).

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use gorilla_settings_mod, only: load_gorilla_inp
    use tetra_grid_settings_mod, only: load_tetra_grid_inp
    use tetra_grid_mod, only: ntetr, verts_rphiz
    use tetra_grid_settings_mod, only: n_field_periods
    use find_tetra_mod, only: find_tetra
    use neutral_flight_mod, only: face_cartesian, face_hybrid, nf_ionised, nf_wall, ho_failed, &
                                  flight_result_t, face_sag_measure, choose_face_mode, fly_neutral, &
                                  lambda_mesh, cart_of_rphiz
    use thermal_force_rng_mod, only: rng_t, rng_seed, rng_uniform
    use marker_species_mod, only: species_set_t, init_species_set, set_marker_charge, ion_to_neutral, &
                                  neutral_to_ion, field_direction

    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
    real(dp), parameter :: r0 = 200.0_dp, a_minor = 50.0_dp, vel = 1.0e6_dp
    integer, parameter :: n_launch = 8, n_stat = 20000
    integer :: n2, nfp, n_fail, mode_auto, mode_expect, i, nb, n_deg
    real(dp) :: sag_tol, ratio, sag, sag_exp, sag_max, per, rmax, dr_max
    real(dp), allocatable :: poly(:,:)
    real(dp) :: xs(3, n_launch), dirs(3, n_launch), t_hyb(n_launch), t_ref, p_hyb(3, n_launch)
    character(len=32) :: arg
    type(flight_result_t) :: res

    n_fail = 0
    call get_command_argument(1, arg); read(arg, *) n2
    call get_command_argument(2, arg); read(arg, *) nfp
    call get_command_argument(3, arg); read(arg, *) sag_tol
    call get_command_argument(4, arg)
    mode_expect = face_hybrid
    if (trim(arg) == 'cartesian') mode_expect = face_cartesian

    call write_inputs(n2, nfp)
    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla()
    per = 2.0_dp*pi/n_field_periods
    print '(a,i0,a,i0,a,i0)', 'mesh: ntetr = ', ntetr, ', n2 = ', n2, ', n_field_periods = ', n_field_periods

    ! 1. chord sag
    call face_sag_measure(ratio, sag, n_deg)
    rmax = maxval(verts_rphiz(1,:))
    sag_exp = rmax*(1.0_dp - cos(0.5_dp*per/n2))
    dr_max = tetra_r_extent()
    sag_max = sag_exp + 0.5_dp*dr_max*per/n2
    print '(a,es12.4,a,es12.4,a,es12.4,a,es12.4)', 'max face sag ', sag, ' cm, chord sag ', sag_exp, &
        ', bound with twist ', sag_max, ', max sag/inradius ', ratio
    print '(a,i0,a)', 'zero-volume tetrahedra skipped (flagged placeholder slots of the grid_kind 5 axis fan): ', n_deg, ''
    if (sag < 0.8_dp*sag_exp .or. sag > sag_max) call fail('face sag')

    ! 2. face_hybrid wall exits vs the polygon
    call boundary_polygon(poly, nb)
    print '(a,i0,a)', 'boundary polygon: ', nb, ' vertices'
    ! launch points (R, phi, Z) and directions (e_R, e_phi, e_Z components)
    xs(:,1) = [220.0_dp, 0.3_dp, 0.0_dp];          dirs(:,1) = [1.0_dp, 0.0_dp, 0.0_dp]
    xs(:,2) = [220.0_dp, 0.3_dp, 10.0_dp];         dirs(:,2) = [0.6_dp, 0.8_dp, 0.0_dp]
    xs(:,3) = [210.0_dp, per - 0.05_dp, -5.0_dp];  dirs(:,3) = [0.3_dp, 0.9_dp, 0.3_dp]   ! period seam
    xs(:,4) = [230.0_dp, 1.0_dp, 20.0_dp];         dirs(:,4) = [-0.5_dp, 0.2_dp, 0.7_dp]
    xs(:,5) = [240.0_dp, 0.5_dp, 10.0_dp];         dirs(:,5) = [-1.0_dp, 0.0_dp, 0.0_dp]  ! across, past the axis ring
    xs(:,6) = [200.0_dp, 0.2_dp, 30.0_dp];         dirs(:,6) = [0.0_dp, 1.0_dp, 0.0_dp]   ! toroidal
    ! through the triangle fan at the magnetic axis (R0 = 200, Z = 0): 0.8 cm from it, and exactly through it
    xs(:,7) = [240.0_dp, 0.5_dp, 0.0_dp];          dirs(:,7) = [-1.0_dp, 0.0_dp, 0.02_dp]
    xs(:,8) = [240.0_dp, 0.5_dp, -0.8_dp];         dirs(:,8) = [-1.0_dp, 0.0_dp, 0.02_dp]
    do i = 1, n_launch
        xs(2,i) = modulo(xs(2,i), per)
        call launch(face_hybrid, i, 0.0_dp, 1.0_dp, res)
        t_hyb(i) = res%t
        p_hyb(:,i) = res%p
        t_ref = polygon_exit_time(i)
        print '(a,i0,a,es14.6,a,es14.6,a,es10.2,a,i0,a)', '  launch ', i, ': hybrid wall time ', res%t, &
            ', polygon ', t_ref, ', |dx| ', abs(res%t - t_ref)*vel, ' cm, ', res%n_cross, ' crossings'
        if (res%status /= nf_wall) call fail('hybrid: no wall hit')
        if (abs(res%t - t_ref)*vel > 1.0e-6_dp) call fail('hybrid wall exit differs from the polygon')
    end do

    ! 3. 'auto' choice, face_cartesian
    call choose_face_mode('auto', sag_tol, mode_auto, ratio, sag, n_deg)
    if (mode_auto /= mode_expect) call fail("'auto' face test")
    print '(a,a)', "'auto' chooses ", merge('cartesian', 'hybrid   ', mode_auto == face_cartesian)
    if (mode_auto == face_cartesian) then
        do i = 1, n_launch
            call launch(face_cartesian, i, 0.0_dp, 1.0_dp, res)
            print '(a,i0,a,es10.2,a)', '  launch ', i, ': cartesian wall exit - hybrid ', &
                norm2(res%p - p_hyb(:,i)), ' cm'
            if (res%status /= nf_wall) call fail('cartesian: no wall hit')
            ! normal offset of the wall <= sag; along a line at incidence angle alpha: sag/cos(alpha)
            if (i /= 6 .and. norm2(res%p - p_hyb(:,i)) > 20.0_dp*sag) call fail('cartesian wall exit too far from hybrid')
        end do
    end if

    ! 4. ionisation statistics
    call statistics(face_hybrid)
    if (mode_auto == face_cartesian) call statistics(face_cartesian)

    ! 5. guiding centre <-> particle
    call round_trip()

    if (n_fail > 0) then
        print '(a,i0,a)', 'test_neutral_flight: ', n_fail, ' check(s) FAILED'
        error stop 1
    end if
    print '(a)', 'test_neutral_flight: all checks passed'

contains

subroutine write_inputs(n2, nfp)

    integer, intent(in) :: n2, nfp
    integer :: iu

    open(newunit=iu, file='tetra_grid.inp', status='replace', action='write')
    write(iu, '(a)') '&TETRA_GRID_NML'
    write(iu, '(a)') ' grid_kind = 5 , n1 = 10 , n3 = 30 ,'
    write(iu, '(a,i0,a,i0,a)') ' n2 = ', n2, ' , boole_n_field_periods = .false. , n_field_periods_manual = ', nfp, ' ,'
    write(iu, '(a)') ' R0_analytic_circ = 200.0 , a_analytic_circ = 50.0 , B0_analytic_circ = 2.0e4 ,'
    write(iu, '(a)') ' q0_analytic_circ = 2.0 , q1_analytic_circ = 0.0 ,'
    write(iu, '(a)') " boole_write_mesh_obj = .false. , filename_mesh_rphiz = 'mesh_rphiz.obj' ,"
    write(iu, '(a)') " filename_mesh_sthetaphi = 'mesh_sthetaphi.obj' ,"
    write(iu, '(a)') '/'
    close(iu)

end subroutine write_inputs

subroutine launch(mode, i, nu, tau_star, res)

    integer, intent(in) :: mode, i
    real(dp), intent(in) :: nu, tau_star
    type(flight_result_t), intent(out) :: res
    real(dp) :: x(3), p0(3), v(3)
    real(dp), allocatable :: nu_all(:)
    integer :: k, iface

    x = xs(:,i)
    call find_tetra(x, vel, vel, k, iface)
    if (k < 1) error stop 'test: launch point not found'
    p0 = cart_of_rphiz(xs(:,i))
    v = velocity(i)
    allocate(nu_all(ntetr))
    nu_all = nu
    call fly_neutral(mode, p0, v, 0.0_dp, k, 1.0_dp, tau_star, nu_all, res)

end subroutine launch

function velocity(i) result(v)

    integer, intent(in) :: i
    real(dp) :: v(3), d(3), c, s

    d = dirs(:,i)/norm2(dirs(:,i))
    c = cos(xs(2,i)); s = sin(xs(2,i))
    v = vel*[d(1)*c - d(2)*s, d(1)*s + d(2)*c, d(3)]

end function velocity

real(dp) function tetra_r_extent() result(dr)

    ! largest R extent of a tetrahedron

    use tetra_grid_mod, only: tetra_grid

    real(dp) :: rr(4)
    integer :: k, iv

    dr = 0.0_dp
    do k = 1, ntetr
        do iv = 1, 4
            rr(iv) = verts_rphiz(1, tetra_grid(k)%ind_knot(iv))
        end do
        dr = max(dr, maxval(rr) - minval(rr))
    end do

end function tetra_r_extent

subroutine boundary_polygon(poly, nb)

    ! (R, Z) of the vertices on r = a, sorted by poloidal angle

    real(dp), allocatable, intent(out) :: poly(:,:)
    integer, intent(out) :: nb
    real(dp), allocatable :: rz(:,:), ang(:)
    real(dp) :: r, tmp(2), ta
    integer :: iv, j, l
    logical :: dup

    allocate(rz(2, size(verts_rphiz, 2)), ang(size(verts_rphiz, 2)))
    nb = 0
    do iv = 1, size(verts_rphiz, 2)
        r = hypot(verts_rphiz(1,iv) - r0, verts_rphiz(3,iv))
        if (r < a_minor*(1.0_dp - 1.0e-6_dp)) cycle
        dup = .false.
        do j = 1, nb
            if (abs(rz(1,j) - verts_rphiz(1,iv)) + abs(rz(2,j) - verts_rphiz(3,iv)) < 1.0e-9_dp) dup = .true.
        end do
        if (dup) cycle
        nb = nb + 1
        rz(:,nb) = [verts_rphiz(1,iv), verts_rphiz(3,iv)]
        ang(nb) = atan2(rz(2,nb), rz(1,nb) - r0)
    end do
    do j = 2, nb
        tmp = rz(:,j); ta = ang(j); l = j - 1
        do while (l >= 1)
            if (ang(l) <= ta) exit
            rz(:,l+1) = rz(:,l); ang(l+1) = ang(l); l = l - 1
        end do
        rz(:,l+1) = tmp; ang(l+1) = ta
    end do
    allocate(poly(2, nb))
    poly = rz(:, 1:nb)

end subroutine boundary_polygon

logical function inside_polygon(rr, zz)

    real(dp), intent(in) :: rr, zz
    integer :: j, jn

    inside_polygon = .true.
    do j = 1, nb
        jn = modulo(j, nb) + 1
        ! counter-clockwise polygon: inside is to the left of each edge
        if ((poly(1,jn) - poly(1,j))*(zz - poly(2,j)) - (poly(2,jn) - poly(2,j))*(rr - poly(1,j)) < 0.0_dp) then
            inside_polygon = .false.
            return
        end if
    end do

end function inside_polygon

function polygon_exit_time(i) result(t)

    integer, intent(in) :: i
    real(dp) :: t, p0(3), v(3), ta, tb, tm, dt
    integer :: it

    p0 = cart_of_rphiz(xs(:,i))
    v = velocity(i)
    dt = 0.01_dp/vel
    ta = 0.0_dp
    do
        tb = ta + dt
        if (.not. inside_line(p0, v, tb)) exit
        ta = tb
    end do
    do it = 1, 100
        tm = 0.5_dp*(ta + tb)
        if (inside_line(p0, v, tm)) then
            ta = tm
        else
            tb = tm
        end if
    end do
    t = 0.5_dp*(ta + tb)

end function polygon_exit_time

logical function inside_line(p0, v, tt)

    real(dp), intent(in) :: p0(3), v(3), tt
    real(dp) :: p(3)

    p = p0 + v*tt
    inside_line = inside_polygon(hypot(p(1), p(2)), p(3))

end function inside_line

subroutine statistics(mode)

    ! constant rate nu with nu t_wall = 1 on launch 1 (radially outward)

    integer, intent(in) :: mode
    type(rng_t) :: rng
    type(flight_result_t) :: r
    real(dp) :: nu, t_wall, p_wall, p_half, f_wall, f_half, x(3), p0(3), v(3), lam(4)
    real(dp), allocatable :: nu_all(:)
    integer :: ip, k, iface, n_wall, n_half, n_bad, n_hof

    t_wall = t_hyb(1)
    nu = 1.0_dp/t_wall
    allocate(nu_all(ntetr))
    nu_all = nu
    x = xs(:,1)
    call find_tetra(x, vel, vel, k, iface)
    p0 = cart_of_rphiz(xs(:,1))
    v = velocity(1)
    n_wall = 0; n_half = 0; n_bad = 0; n_hof = 0
    do ip = 1, n_stat
        call rng_seed(rng, 7_int64, int(ip, int64))
        call fly_neutral(mode, p0, v, 0.0_dp, k, 1.0_dp, -log(1.0_dp - rng_uniform(rng)), nu_all, r)
        if (r%status == nf_wall) n_wall = n_wall + 1
        if (r%status == nf_ionised) then
            if (r%t < 0.5_dp*t_wall) n_half = n_half + 1
            if (r%handover == ho_failed) then
                n_hof = n_hof + 1
            else
                lam = lambda_mesh(r%ind_tetr, r%x)
                if (any(lam < -1.0e-8_dp)) n_bad = n_bad + 1
            end if
        end if
    end do
    p_wall = exp(-1.0_dp); p_half = 1.0_dp - exp(-0.5_dp)
    f_wall = real(n_wall, dp)/n_stat; f_half = real(n_half, dp)/n_stat
    print '(a,a,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,i0,a,i0)', ' statistics (', trim(merge('cartesian', 'hybrid   ', &
        mode == face_cartesian)), '): wall ', f_wall, ' (exact ', p_wall, '), ionised before t_wall/2 ', f_half, &
        ' (exact ', p_half, '); outside returned tetrahedron ', n_bad, ', handover failed ', n_hof
    if (abs(f_wall - p_wall) > 4.5_dp*sqrt(p_wall*(1.0_dp - p_wall)/n_stat)) call fail('wall fraction')
    if (abs(f_half - p_half) > 4.5_dp*sqrt(p_half*(1.0_dp - p_half)/n_stat)) call fail('ionisation time distribution')
    if (n_bad > 0 .or. n_hof > 0) call fail('ionisation point not in the returned tetrahedron')

end subroutine statistics

subroutine round_trip()

    use constants, only: echarge, clight

    real(dp), parameter :: amu = 1.66053906660e-24_dp
    type(species_set_t) :: ss
    real(dp) :: rate0(0:4), xg(3), x(3), vcart(3), jump(3), jump2(3), vpar, vperp, vpar2, vperp2, b(3), bmod, &
                rho, err_x, err_v, err_perp, err_len
    integer :: k0, k, iface, i

    rate0 = 0.0_dp
    call init_species_set(ss, 'Be', [2.01410177812_dp*amu], [1.0_dp], [1.0e13_dp], [20.0_dp], 20.0_dp, &
                          rate0, rate0, .true.)
    call set_marker_charge(ss, 1)
    xg = [230.0_dp, 0.4_dp, 5.0_dp]
    call find_tetra(xg, vel, vel, k0, iface)
    vpar = 1.0e6_dp; vperp = 2.0e6_dp
    call field_direction(xg, k0, b, bmod)
    rho = vperp*ss%mass_g*clight/(echarge*bmod)
    err_x = 0.0_dp; err_v = 0.0_dp; err_perp = 0.0_dp; err_len = 0.0_dp
    do i = 0, 7
        x = xg; k = k0; iface = 0
        call ion_to_neutral(ss, 1, (i + 0.3_dp)/8.0_dp, x, k, iface, vpar, vperp, vcart, jump)
        err_len = max(err_len, abs(norm2(jump) - rho)/rho)
        err_perp = max(err_perp, abs(dot_product(jump, b))/rho)
        err_v = max(err_v, abs(norm2(vcart) - hypot(vpar, vperp))/hypot(vpar, vperp))
        call neutral_to_ion(ss, x, k, iface, vcart, vpar2, vperp2, jump2)
        err_x = max(err_x, norm2(cart_of_rphiz(x) - cart_of_rphiz(xg))/rho)
        err_v = max(err_v, abs(vpar2 - vpar)/vperp, abs(vperp2 - vperp)/vperp)
    end do
    print '(a,es10.3,a,4es10.2)', ' round trip (Be1+, Larmor radius ', rho, ' cm): |X''-X|/rho, |jump|-rho, ' // &
        'jump.b/rho, velocity error: ', err_x, err_len, err_perp, err_v
    if (err_len > 1.0e-10_dp .or. err_perp > 1.0e-10_dp) call fail('ion_to_neutral: Larmor displacement')
    if (err_x > 1.0e-2_dp .or. err_v > 1.0e-2_dp) call fail('round trip ion -> neutral -> ion')

end subroutine round_trip

subroutine fail(msg)

    character(len=*), intent(in) :: msg

    print '(a,a)', 'FAIL: ', msg
    n_fail = n_fail + 1

end subroutine fail

end program test_neutral_flight
