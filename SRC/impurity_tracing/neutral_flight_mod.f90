module neutral_flight_mod

! Straight-line flight of neutrals through the GORILLA tetrahedral mesh (coord_system = 1, R phi Z).
!
! The neutral moves on a straight line in Cartesian lab coordinates, p(t) = p0 + v (t - t0). For each
! tetrahedron the line is rotated into the field period of that tetrahedron (multiples of
! 2 pi/n_field_periods), so periodic boundaries need no special treatment. The exit face is found with
! one of two face tests (barycentric coordinates lambda_i: 0 on face i, 1 at the opposite vertex):
!   face_hybrid     faces of the GORILLA tetrahedron, planar in (R, phi, Z), evaluated along the
!                   Cartesian line; exit time by bracketing and safeguarded Newton. Exact for any
!                   number of toroidal slices; the neutral is always in the tetrahedron the guiding-
!                   centre pusher uses, so the handover at ionisation needs no search.
!   face_cartesian  planar faces through the Cartesian images of the vertices; exit time is linear.
!                   The toroidal arcs become chords (sag R (1 - cos(dphi/2)) for a slice dphi), so it
!                   needs fine toroidal slices; at ionisation the cylindrical tetrahedron is searched
!                   (same tetrahedron, its face neighbours, find_tetra).
! face_sag_measure gives the largest chord sag relative to the tetrahedron inradius over the mesh;
! choose_face_mode('auto') uses face_cartesian if it is below a tolerance, else face_hybrid.
!
! Ionisation: track-length estimator. The caller draws tau_star = -ln(1 - u); the neutral ionises
! where the integral of nu_tetr(ind_tetr) dt along the path reaches tau_star (rate constant per
! tetrahedron, exact). Rotating frames are not implemented yet (the line is then no longer straight in
! the mesh frame; face_hybrid extends to it by adding -Omega t to phi).
! CGS units.

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
    real(dp), parameter :: lam_tol = 1.0e-9_dp      ! inside test, barycentric
    integer, parameter :: max_cross = 1000000

    integer, parameter, public :: face_cartesian = 1, face_hybrid = 2
    ! flight status
    integer, parameter, public :: nf_ionised = 1, nf_wall = 2, nf_time = 3, nf_error = 4
    ! how the cylindrical tetrahedron was found at the end of a face_cartesian flight
    integer, parameter, public :: ho_direct = 0, ho_same = 1, ho_neighbour = 2, ho_find = 3, ho_failed = 4

    type, public :: flight_result_t
        integer :: status = nf_error
        real(dp) :: t = 0.0_dp          ! time at the end of the flight [s]
        real(dp) :: p(3) = 0.0_dp       ! Cartesian lab position at the end [cm]
        real(dp) :: x(3) = 0.0_dp       ! (R, phi, Z) at the end, phi in the range of ind_tetr
        integer :: ind_tetr = -1        ! tetrahedron at the end (the last one inside for nf_wall)
        integer :: n_cross = 0          ! tetrahedron crossings
        integer :: handover = ho_direct
    end type flight_result_t

    public :: face_sag_measure, choose_face_mode, fly_neutral, exit_tetrahedron, locate_point, &
              lambda_mesh, cart_of_rphiz

contains

!-------------------------------------------------------------------------------------------------
! geometry helpers

pure function period_length() result(per)

    use tetra_grid_settings_mod, only: n_field_periods

    real(dp) :: per

    per = 2.0_dp*pi/real(n_field_periods, dp)

end function period_length

pure function cart_of_rphiz(x) result(p)

    real(dp), intent(in) :: x(3)
    real(dp) :: p(3)

    p = [x(1)*cos(x(2)), x(1)*sin(x(2)), x(3)]

end function cart_of_rphiz

subroutine tetra_vertices(k, xv)

    ! (R, phi, Z) of the 4 vertices of tetrahedron k, phi unwrapped to the range of the tetrahedron
    ! (next to the first vertex, tetra_physics%x1), in GORILLA's vertex order

    use tetra_grid_mod, only: tetra_grid, verts_rphiz
    use tetra_physics_mod, only: tetra_physics

    integer, intent(in) :: k
    real(dp), intent(out) :: xv(3,4)
    real(dp) :: per, phic
    integer :: iv

    per = period_length()
    phic = tetra_physics(k)%x1(2)
    do iv = 1, 4
        xv(:,iv) = verts_rphiz(:, tetra_grid(k)%ind_knot(iv))
        xv(2,iv) = xv(2,iv) + per*anint((phic - xv(2,iv))/per)
    end do

end subroutine tetra_vertices

subroutine to_tetra_frame(k, p, v, pk, vk)

    ! Rotate the Cartesian point p and direction v about the Z axis by a multiple of the field period
    ! so that the toroidal angle of p is next to the range of tetrahedron k

    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_settings_mod, only: n_field_periods

    integer, intent(in) :: k
    real(dp), intent(in) :: p(3), v(3)
    real(dp), intent(out) :: pk(3), vk(3)
    real(dp) :: per, alpha, c, s
    integer :: m

    per = period_length()
    m = nint((tetra_physics(k)%x1(2) - atan2(p(2), p(1)))/per)
    m = modulo(m, n_field_periods)
    if (m == 0) then
        pk = p; vk = v
        return
    end if
    alpha = m*per
    c = cos(alpha); s = sin(alpha)
    pk = [c*p(1) - s*p(2), s*p(1) + c*p(2), p(3)]
    vk = [c*v(1) - s*v(2), s*v(1) + c*v(2), v(3)]

end subroutine to_tetra_frame

pure function rphiz_near(p, phic) result(x)

    ! (R, phi, Z) of the Cartesian point p with phi unwrapped next to phic (same frame, no rotation)

    real(dp), intent(in) :: p(3), phic
    real(dp) :: x(3)

    x(1) = sqrt(p(1)**2 + p(2)**2)
    x(2) = atan2(p(2), p(1))
    x(2) = x(2) + 2.0_dp*pi*anint((phic - x(2))/(2.0_dp*pi))
    x(3) = p(3)

end function rphiz_near

function lambda_mesh(k, x) result(lam)

    ! Barycentric coordinates of the (R, phi, Z) point x in tetrahedron k (GORILLA faces):
    ! lam(i) = normal distance to face i / normal distance of vertex i

    use tetra_physics_mod, only: tetra_physics

    integer, intent(in) :: k
    real(dp), intent(in) :: x(3)
    real(dp) :: lam(4), z(3)

    z = x - tetra_physics(k)%x1
    lam = matmul(z, tetra_physics(k)%anorm)
    lam(1) = lam(1) + tetra_physics(k)%dist_ref
    lam = lam/tetra_physics(k)%dist_ref_vec

end function lambda_mesh

subroutine cartesian_planes(k, nrm, off)

    ! Planar faces through the Cartesian images of the vertices of k: lam(i) = nrm(:,i).p + off(i)

    integer, intent(in) :: k
    real(dp), intent(out) :: nrm(3,4), off(4)
    real(dp) :: xv(3,4), c(3,4), a(3), e1(3), e2(3), n(3), den
    integer :: i, j, l, idx(3)

    call tetra_vertices(k, xv)
    do i = 1, 4
        c(:,i) = cart_of_rphiz(xv(:,i))
    end do
    do i = 1, 4
        l = 0
        do j = 1, 4
            if (j == i) cycle
            l = l + 1
            idx(l) = j
        end do
        a = c(:,idx(1)); e1 = c(:,idx(2)) - a; e2 = c(:,idx(3)) - a
        n = [e1(2)*e2(3) - e1(3)*e2(2), e1(3)*e2(1) - e1(1)*e2(3), e1(1)*e2(2) - e1(2)*e2(1)]
        den = dot_product(n, c(:,i) - a)
        nrm(:,i) = n/den
        off(i) = -dot_product(n, a)/den
    end do

end subroutine cartesian_planes

!-------------------------------------------------------------------------------------------------
! mode selection

subroutine face_sag_measure(max_ratio, max_sag, n_degenerate)

    ! Largest distance of the Cartesian images of face edge midpoints and face centroids (taken in
    ! (R, phi, Z)) from the planar Cartesian face through the face vertices: max_sag [cm], and
    ! max_ratio = largest such distance / inradius of the Cartesian tetrahedron. Tetrahedra of zero
    ! volume (e.g. collapsed onto the magnetic axis) are skipped and counted in n_degenerate.

    use tetra_grid_mod, only: ntetr

    real(dp), intent(out) :: max_ratio, max_sag
    integer, intent(out) :: n_degenerate
    real(dp) :: xv(3,4), c(3,4), a(3), e1(3), e2(3), n(3), area(4), vol, rin, sag, xm(3), diam
    integer :: k, i, j, l, m, idx(3)

    max_ratio = 0.0_dp; max_sag = 0.0_dp; n_degenerate = 0
    !$omp parallel do default(shared) private(k, i, j, l, m, idx, xv, c, a, e1, e2, n, area, vol, rin, sag, xm, diam) &
    !$omp& reduction(max:max_ratio, max_sag) reduction(+:n_degenerate)
    do k = 1, ntetr
        call tetra_vertices(k, xv)
        do i = 1, 4
            c(:,i) = cart_of_rphiz(xv(:,i))
        end do
        e1 = c(:,2) - c(:,1); e2 = c(:,3) - c(:,1)
        n = [e1(2)*e2(3) - e1(3)*e2(2), e1(3)*e2(1) - e1(1)*e2(3), e1(1)*e2(2) - e1(2)*e2(1)]
        vol = abs(dot_product(n, c(:,4) - c(:,1)))/6.0_dp
        diam = 0.0_dp
        do i = 1, 3
            do j = i + 1, 4
                diam = max(diam, norm2(c(:,i) - c(:,j)))
            end do
        end do
        if (vol <= 1.0e-10_dp*diam**3) then
            n_degenerate = n_degenerate + 1
            cycle
        end if
        sag = 0.0_dp
        do i = 1, 4
            l = 0
            do j = 1, 4
                if (j == i) cycle
                l = l + 1
                idx(l) = j
            end do
            a = c(:,idx(1)); e1 = c(:,idx(2)) - a; e2 = c(:,idx(3)) - a
            n = [e1(2)*e2(3) - e1(3)*e2(2), e1(3)*e2(1) - e1(1)*e2(3), e1(1)*e2(2) - e1(2)*e2(1)]
            area(i) = 0.5_dp*norm2(n)
            n = n/norm2(n)
            do m = 1, 4
                select case (m)
                case (1); xm = 0.5_dp*(xv(:,idx(1)) + xv(:,idx(2)))
                case (2); xm = 0.5_dp*(xv(:,idx(2)) + xv(:,idx(3)))
                case (3); xm = 0.5_dp*(xv(:,idx(3)) + xv(:,idx(1)))
                case (4); xm = (xv(:,idx(1)) + xv(:,idx(2)) + xv(:,idx(3)))/3.0_dp
                end select
                sag = max(sag, abs(dot_product(n, cart_of_rphiz(xm) - a)))
            end do
        end do
        rin = 3.0_dp*vol/sum(area)
        max_sag = max(max_sag, sag)
        if (rin > 0.0_dp) max_ratio = max(max_ratio, sag/rin)
    end do
    !$omp end parallel do

end subroutine face_sag_measure

subroutine choose_face_mode(mode_name, sag_tol, mode, max_ratio, max_sag, n_degenerate)

    ! mode_name: 'auto', 'cartesian' or 'hybrid'. 'auto': face_cartesian if max_ratio <= sag_tol.

    use tetra_physics_mod, only: coord_system

    character(len=*), intent(in) :: mode_name
    real(dp), intent(in) :: sag_tol
    integer, intent(out) :: mode
    real(dp), intent(out) :: max_ratio, max_sag
    integer, intent(out) :: n_degenerate

    if (coord_system /= 1) error stop 'neutral flight: needs coord_system = 1 (R, phi, Z)'
    call face_sag_measure(max_ratio, max_sag, n_degenerate)
    select case (trim(adjustl(mode_name)))
    case ('auto')
        mode = face_hybrid
        if (max_ratio <= sag_tol) mode = face_cartesian
    case ('cartesian')
        mode = face_cartesian
        if (max_ratio > sag_tol) print '(a,es10.3,a,es10.3,a)', ' neutral flight: WARNING face_cartesian with chord sag ', &
            max_ratio, ' of the inradius (tolerance ', sag_tol, '): use more toroidal slices or hybrid'
    case ('hybrid')
        mode = face_hybrid
    case default
        error stop "neutral flight: face test must be 'auto', 'cartesian' or 'hybrid'"
    end select

end subroutine choose_face_mode

!-------------------------------------------------------------------------------------------------
! exit of one tetrahedron

subroutine eval_hybrid(k, pk, vk, dt, lam, dlam)

    ! lam and d lam/dt of the GORILLA faces of k at the point pk + vk dt (tetrahedron frame)

    use tetra_physics_mod, only: tetra_physics

    integer, intent(in) :: k
    real(dp), intent(in) :: pk(3), vk(3), dt
    real(dp), intent(out) :: lam(4), dlam(4)
    real(dp) :: p(3), x(3), dx(3), r2

    p = pk + vk*dt
    x = rphiz_near(p, tetra_physics(k)%x1(2))
    r2 = p(1)**2 + p(2)**2
    dx = [(p(1)*vk(1) + p(2)*vk(2))/x(1), (p(1)*vk(2) - p(2)*vk(1))/r2, vk(3)]
    lam = lambda_mesh(k, x)
    dlam = matmul(dx, tetra_physics(k)%anorm)/tetra_physics(k)%dist_ref_vec

end subroutine eval_hybrid

subroutine exit_tetrahedron(mode, k, iface_in, p0, v, dt_exit, iface_out)

    ! Time dt_exit >= 0 after which the line p0 + v dt (lab frame) leaves tetrahedron k, and the exit
    ! face. iface_in: entry face (excluded from the first bracket), 0 if the line starts inside.

    use tetra_physics_mod, only: tetra_physics

    integer, intent(in) :: mode, k, iface_in
    real(dp), intent(in) :: p0(3), v(3)
    real(dp), intent(out) :: dt_exit
    integer, intent(out) :: iface_out
    real(dp) :: pk(3), vk(3), nrm(3,4), off(4), lam(4), dlam(4), lam_b(4), dlam_b(4), ta, tb, tr, step, diam
    real(dp) :: xv(3,4), c(3,4)
    integer :: i, j, iter
    logical :: first

    call to_tetra_frame(k, p0, v, pk, vk)
    dt_exit = huge(1.0_dp)
    iface_out = 0
    if (any(tetra_physics(k)%dist_ref_vec <= 0.0_dp)) then
        print '(a,i0,a)', ' neutral flight: tetrahedron ', k, ' has zero volume'
        error stop 'neutral flight: degenerate tetrahedron on the path'
    end if

    if (mode == face_cartesian) then
        call cartesian_planes(k, nrm, off)
        do i = 1, 4
            if (i == iface_in) cycle
            lam(i) = dot_product(nrm(:,i), pk) + off(i)
            dlam(i) = dot_product(nrm(:,i), vk)
            if (dlam(i) < 0.0_dp) then
                tr = max(0.0_dp, -lam(i)/dlam(i))
                if (tr < dt_exit) then
                    dt_exit = tr; iface_out = i
                end if
            end if
        end do
        if (iface_out == 0) error stop 'neutral flight: no exit face (face_cartesian)'
        return
    end if

    ! face_hybrid: march with the linear estimate of the next crossing until a face changes sign,
    ! then find the first root in the bracket
    call tetra_vertices(k, xv)
    do i = 1, 4
        c(:,i) = cart_of_rphiz(xv(:,i))
    end do
    diam = 0.0_dp
    do i = 1, 3
        do j = i + 1, 4
            diam = max(diam, norm2(c(:,i) - c(:,j)))
        end do
    end do
    ta = 0.0_dp
    call eval_hybrid(k, pk, vk, ta, lam, dlam)
    first = .true.
    do iter = 1, 200
        step = huge(1.0_dp)
        do i = 1, 4
            if (first .and. i == iface_in) cycle
            if (dlam(i) < 0.0_dp) step = min(step, max(lam(i), 0.0_dp)/(-dlam(i)))
        end do
        step = min(1.02_dp*step + 1.0e-6_dp*diam/norm2(vk), 1.1_dp*diam/norm2(vk))
        tb = ta + step
        call eval_hybrid(k, pk, vk, tb, lam_b, dlam_b)
        do i = 1, 4
            if (first .and. i == iface_in) cycle
            if (lam_b(i) < 0.0_dp) then
                tr = root_in_bracket(k, i, pk, vk, ta, tb)
                if (tr < dt_exit) then
                    dt_exit = tr; iface_out = i
                end if
            end if
        end do
        if (iface_out /= 0) return
        ta = tb; lam = lam_b; dlam = dlam_b
        first = .false.
    end do
    error stop 'neutral flight: no exit face found (face_hybrid)'

end subroutine exit_tetrahedron

function root_in_bracket(k, i, pk, vk, ta_in, tb_in) result(tr)

    ! Root of lam_i(t) in [ta, tb] with lam_i(ta) >= 0 > lam_i(tb): Newton, bisection as fallback

    integer, intent(in) :: k, i
    real(dp), intent(in) :: pk(3), vk(3), ta_in, tb_in
    real(dp) :: tr, ta, tb, lam(4), dlam(4), tn
    integer :: iter

    ta = ta_in; tb = tb_in
    tr = 0.5_dp*(ta + tb)
    do iter = 1, 100
        call eval_hybrid(k, pk, vk, tr, lam, dlam)
        if (lam(i) >= 0.0_dp) then
            ta = tr
        else
            tb = tr
        end if
        if (abs(lam(i)) < 1.0e-14_dp .or. tb - ta <= 1.0e-15_dp*max(abs(tb), 1.0e-30_dp)) exit
        tn = tr - lam(i)/dlam(i)
        if (dlam(i) == 0.0_dp .or. tn <= ta .or. tn >= tb) tn = 0.5_dp*(ta + tb)
        tr = tn
    end do

end function root_in_bracket

!-------------------------------------------------------------------------------------------------

subroutine locate_point(k, p, vnorm, kout, x, how)

    ! Cylindrical (GORILLA) tetrahedron of the Cartesian lab point p, starting the search at k:
    ! k itself, its face neighbours, then find_tetra. x: (R, phi, Z) of p in the range of kout.

    use tetra_grid_mod, only: tetra_grid
    use tetra_physics_mod, only: tetra_physics
    use tetra_grid_settings_mod, only: n_field_periods
    use find_tetra_mod, only: find_tetra

    integer, intent(in) :: k
    real(dp), intent(in) :: p(3), vnorm
    integer, intent(out) :: kout, how
    real(dp), intent(out) :: x(3)
    real(dp) :: xs(3), per
    integer :: j, kn, iface

    if (point_in(k, p, x)) then
        kout = k; how = ho_same
        return
    end if
    do j = 1, 4
        kn = tetra_grid(k)%neighbour_tetr(j)
        if (kn < 1) cycle
        if (point_in(kn, p, x)) then
            kout = kn; how = ho_neighbour
            return
        end if
    end do
    per = 2.0_dp*pi/real(n_field_periods, dp)
    xs(1) = sqrt(p(1)**2 + p(2)**2)
    xs(2) = modulo(atan2(p(2), p(1)), per)
    xs(3) = p(3)
    call find_tetra(xs, vnorm, vnorm, kout, iface)
    x = xs
    how = ho_find
    if (kout < 1) how = ho_failed

contains

    logical function point_in(kk, pp, xx)
        integer, intent(in) :: kk
        real(dp), intent(in) :: pp(3)
        real(dp), intent(out) :: xx(3)
        real(dp) :: pk(3), dummy(3)
        call to_tetra_frame(kk, pp, [0.0_dp, 0.0_dp, 0.0_dp], pk, dummy)
        xx = rphiz_near(pk, tetra_physics(kk)%x1(2))
        point_in = all(lambda_mesh(kk, xx) >= -lam_tol)
    end function point_in

end subroutine locate_point

subroutine fly_neutral(mode, p0, v, t0, ind_tetr, t_max, tau_star, nu_tetr, res, tri_of_tetr, resid)

    ! Fly a neutral from the Cartesian lab position p0 at time t0 with velocity v [cm/s], starting in
    ! tetrahedron ind_tetr (GORILLA's, containing p0), until it ionises (accumulated optical depth
    ! tau_star with the rates nu_tetr(ntetr) [1/s]), leaves the mesh, or reaches t_max.
    ! Optional: resid(tri_of_tetr(k)) += time spent in tetrahedron k.

    integer, intent(in) :: mode, ind_tetr
    real(dp), intent(in) :: p0(3), v(3), t0, t_max, tau_star, nu_tetr(:)
    type(flight_result_t), intent(out) :: res
    integer, intent(in), optional :: tri_of_tetr(:)
    real(dp), intent(inout), optional :: resid(:)
    real(dp) :: t, tau, dt_exit, dt_seg, nu, pk(3), vk(3), vnorm
    integer :: k, iface_in, iface_out, kn, n, how
    logical :: boole_res

    boole_res = present(tri_of_tetr) .and. present(resid)
    vnorm = norm2(v)
    k = ind_tetr
    t = t0
    tau = 0.0_dp
    iface_in = 0
    n = 0
    if (mode == face_cartesian) then
        ! starting tetrahedron of the planar Cartesian mesh: k or a face neighbour
        call locate_cartesian(k)
    end if

    do
        call exit_tetrahedron(mode, k, iface_in, p0 + v*(t - t0), v, dt_exit, iface_out)
        nu = nu_tetr(k)
        dt_seg = dt_exit
        res%status = 0
        if (nu > 0.0_dp) then
            if (tau + nu*dt_exit >= tau_star) then
                dt_seg = max(0.0_dp, (tau_star - tau)/nu)
                res%status = nf_ionised
            end if
        end if
        if (t + dt_seg >= t_max) then
            dt_seg = max(0.0_dp, t_max - t)
            res%status = nf_time
        end if
        if (boole_res) resid(tri_of_tetr(k)) = resid(tri_of_tetr(k)) + dt_seg
        tau = tau + nu*dt_seg
        t = t + dt_seg
        if (res%status /= 0) exit
        kn = tetra_grid_neighbour(k, iface_out)
        if (kn < 1) then
            res%status = nf_wall
            exit
        end if
        iface_in = tetra_grid_neighbour_face(k, iface_out)
        k = kn
        n = n + 1
        if (n > max_cross) then
            res%status = nf_error
            exit
        end if
    end do

    res%t = t
    res%p = p0 + v*(t - t0)
    res%n_cross = n
    if (mode == face_hybrid .or. res%status == nf_wall) then
        call to_tetra_frame(k, res%p, v, pk, vk)
        res%x = rphiz_near(pk, tetra_x1_phi(k))
        res%ind_tetr = k
        res%handover = ho_direct
    else
        call locate_point(k, res%p, vnorm, res%ind_tetr, res%x, how)
        res%handover = how
    end if

contains

    subroutine locate_cartesian(kk)
        integer, intent(inout) :: kk
        real(dp) :: nrm(3,4), off(4), q(3), dum(3)
        integer :: j, kc
        call cartesian_planes(kk, nrm, off)
        call to_tetra_frame(kk, p0, v, q, dum)
        if (all(matmul(q, nrm) + off >= -lam_tol)) return
        do j = 1, 4
            kc = tetra_grid_neighbour(kk, j)
            if (kc < 1) cycle
            call cartesian_planes(kc, nrm, off)
            call to_tetra_frame(kc, p0, v, q, dum)
            if (all(matmul(q, nrm) + off >= -lam_tol)) then
                kk = kc
                return
            end if
        end do
        ! not found: keep kk (start very close to a face)
    end subroutine locate_cartesian

end subroutine fly_neutral

integer function tetra_grid_neighbour(k, i)
    use tetra_grid_mod, only: tetra_grid
    integer, intent(in) :: k, i
    tetra_grid_neighbour = tetra_grid(k)%neighbour_tetr(i)
end function tetra_grid_neighbour

integer function tetra_grid_neighbour_face(k, i)
    use tetra_grid_mod, only: tetra_grid
    integer, intent(in) :: k, i
    tetra_grid_neighbour_face = tetra_grid(k)%neighbour_face(i)
end function tetra_grid_neighbour_face

real(dp) function tetra_x1_phi(k)
    use tetra_physics_mod, only: tetra_physics
    integer, intent(in) :: k
    tetra_x1_phi = tetra_physics(k)%x1(2)
end function tetra_x1_phi

end module neutral_flight_mod
