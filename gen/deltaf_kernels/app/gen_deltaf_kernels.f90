program gen_deltaf_kernels
    ! FortSym generator for the pointwise kernels of the linear delta-f NTV
    ! estimator. Output (committed): SRC/generated/. Regenerate with
    !     cd gen/deltaf_kernels && fpm run gen_deltaf_kernels -- ../../SRC/generated
    ! and check with gen/deltaf_kernels/check_generated.sh.
    !
    ! Gaussian units. mu = m v_perp**2/(2 B) = -particle_mass*perpinv.
    ! Model: A_phi, covariant B_k axisymmetric; the perturbation enters via
    ! |B| and Phi only (Boozer chartmap, grid_kind = 6), no delta A_par.
    !   p_phi = e A_phi/c + m v_par h_phi,
    !   dp_phi/dt = -(mu + m v_par**2/B) dB/dphi - e dPhi/dphi = -dH1/dphi  [P1, P2]
    ! (mu gB(3) + e gPhi(3) in GORILLA is dH/dphi at fixed (x, v_par); the
    ! m v_par**2/B term comes from dh_k/dphi.) Background: canonical Maxwellian
    ! f0 = F(p_phi, H0) with profiles at psi* = c p_phi/e, so {f0, H0} = 0 [P3].
    ! Along unperturbed orbits -{f0,H1} = F_p dH1/dphi + F_H dH1/dt|0 [P4]; the
    ! adiabatic second term gives no torque, so the weight equation is
    !   dw/dt = d ln f0/dp_phi|_H0 * dH1/dphi,
    !   d ln f0/dp_phi = (c/e) [n'/n + e Phi0'/T + ((H0 - e Phi0)/T - 3/2) T'/T]  [P5],
    ! collisions omitted, and the torque estimator is sum_i w_i <-dH1/dphi>_i.
    ! Circular chart (grid_kind = 5): R = R0 + r cos(theta), q = q0 + q1 (r/a)**2,
    ! psi_pol = B0 a**2/(2 q1) ln(q/q0), B = (B0/R) sqrt(R0**2 + r**2/q**2) [C1, C2].
    ! Harmonic: delta B, delta Phi = (b_mn, Phi_mn)(s) cos(m theta - n phi).
    !
    ! ------------------------------------------------------------------------
    ! Collisional background drive (opt-in, canonical_ntv_mod ntv_fax_file)
    ! ------------------------------------------------------------------------
    ! The collisional axisymmetric background is f_M + f_AX, not the canonical
    ! Maxwellian f0 alone (Kernbichler et al. PPCF 58 (2016) 104001 Eq. (27),
    ! Kasilov et al. PoP 21 (2014) 092506 Eq. (65)). Split f_AX = G + h with
    ! the canonical shift G = f0 - f_M = m v_par h_phi df0/dp_phi (to first
    ! order in rho*), and write the collisional remainder as
    !   h/f0 = (d ln f0/dp_phi) m v h_phi chi(position, lambda),
    ! i.e. chi = lambda h/G with lambda = v_par/v, so chi = lambda for h = G.
    ! chi is imported from the axisymmetric NEO-2 solution (thin orbit, local).
    ! The canonical weight source -{f0, H1}/f0 already contains G; the missing
    ! part is -{h, H1}/f0. In (x, v, lambda) h is axisymmetric, v is conserved
    ! (magnetostatic H1 = mu dB) and the perturbed mirror force is [F1]
    !   d lambda/dt = -(mu/(m v)) b.grad B = -(v (1-lambda**2)/2) b.grad ln B,
    !   delta(d lambda/dt) = -(v (1-lambda**2)/2) b.grad(dB/B).
    ! The perturbed radial drift acting on h is O(rho*) smaller and dropped
    ! (thin orbit, as in NEO-2). For an m = 0 harmonic of the axisymmetric
    ! field b.grad(dB/B) = b^phi d(dB)/dphi / B with b^phi = h_phi/R**2. Then
    !   dw/dt|_fax = -delta(d lambda/dt) d(h/f0)/d lambda
    !              = (d ln f0/dp_phi) m v**2 (h_phi/R)**2 (1-lambda**2)/2
    !                * d chi/d lambda * (d dB/dphi)/B,
    ! with d chi/d lambda at fixed position from a table [F2] (finite at the
    ! turning point lambda = 0, unlike d chi/d eta). This is the mirror term of
    ! Kernbichler Eq. (28) on h, -sigma eta d/d eta = (1-lambda**2)/(2|lambda|)
    ! d/d lambda (NEO-2 f0_ttmp). Check [F3]: for h = G (chi = lambda)
    ! the kernel returns (d ln f0/dp_phi) mu (h_phi/R)**2 d(dB)/dphi, the mu
    ! part of the canonical source (d ln f0/dp_phi) dH1/dphi: the canonical
    ! drive is the radial drift on f_M plus exactly this mirror term on G.
    ! Position label (convention free): vartheta = branch arccos(1 - 2 (B -
    ! Bmin)/(Bmax - Bmin)), branch = sign of dB along b (motion with lambda > 0).
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_string, only: str_t, str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, func, pi_expr, operator(+), &
        operator(-), operator(*), operator(/), operator(**), sin, cos, log, &
        sqrt, acos
    use fortsym_diff, only: diff, partial_derivative
    use fortsym_subs, only: subs_many
    use fortsym_print, only: print_expr
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    implicit none

    character(*), parameter :: FORTSYM_REVISION = &
        "fortsym@8a418e2efe40dd9a01231a0900f71caad7858360"
    character(*), parameter :: GENERATOR = &
        "gen/deltaf_kernels/app/gen_deltaf_kernels.f90"
    character(*), parameter :: REGENERATE = &
        "cd gen/deltaf_kernels && fpm run gen_deltaf_kernels -- ../../SRC/generated"

    type(arena_t), target :: arena
    type(native_engine_t) :: native
    type(symengine_engine_t) :: symeng
    type(suite_t) :: suite
    character(len=4096) :: output_dir

    ! Physical symbols.
    type(expr_t) :: mass, charge, c_light, mu, vpar
    ! Harmonic.
    type(expr_t) :: m_pol, n_tor, theta, phi, s, alpha
    type(expr_t) :: bmn, dbmn_ds, phimn, dphimn_ds
    type(expr_t) :: b_of_s, phi_of_s, db1, dphi1
    ! Circular chart.
    type(expr_t) :: r0, a_minor, b0, q0, q1, r, q, r_major, psi_pol, bmod_c
    type(expr_t) :: b_theta_c, b_phi_c
    type(expr_t) :: bmod
    ! Canonical Maxwellian.
    type(expr_t) :: h0, p_phi, dens, ddens, temp, dtemp, phi0, dphi0
    type(expr_t) :: n_f, t_f, phi0_f, lnf0, dlnf0_closed, dlnf0
    type(expr_t) :: db1_dphi, dphi1_dphi, dh1_dphi
    ! Canonical pairs for the bracket identities.
    type(expr_t) :: th, pth, ph, pph
    ! Collisional background drive.
    type(expr_t) :: bmin, bmax, branch, dchi, h_phi, r_maj, kappa_f

    if (command_argument_count() /= 1) then
        write (output_unit, "(a)") "usage: gen_deltaf_kernels OUTPUT_DIR"
        error stop 2
    end if
    call get_command_argument(1, output_dir)

    call arena%init()
    native = make_native_engine(arena)
    symeng = make_symengine_engine(arena)

    mass = sym(arena, "mass")
    charge = sym(arena, "charge")
    c_light = sym(arena, "c_light")
    mu = sym(arena, "mu")
    vpar = sym(arena, "vpar")
    m_pol = sym(arena, "m_pol")
    n_tor = sym(arena, "n_tor")
    theta = sym(arena, "theta")
    phi = sym(arena, "phi")
    s = sym(arena, "s")
    bmn = sym(arena, "bmn")
    dbmn_ds = sym(arena, "dbmn_ds")
    phimn = sym(arena, "phimn")
    dphimn_ds = sym(arena, "dphimn_ds")

    call suite_begin(suite, "GORILLA_APPLETS linear delta-f kernels")

    ! ---------------------------------------------------------------- harmonic
    alpha = m_pol*theta - n_tor*phi
    b_of_s = func("b_mn", [s])
    phi_of_s = func("Phi_mn", [s])
    db1 = b_of_s*cos(alpha)
    dphi1 = phi_of_s*cos(alpha)
    call emit_harmonic()
    call emit_complex_harmonic()

    ! ---------------------------------------------------------------- circular
    r0 = sym(arena, "R0")
    a_minor = sym(arena, "a_minor")
    b0 = sym(arena, "B0")
    q0 = sym(arena, "q0")
    q1 = sym(arena, "q1")
    r = sym(arena, "r")
    q = q0 + q1*(r/a_minor)**2
    r_major = r0 + r*cos(theta)
    psi_pol = b0*a_minor**2/(2*q1)*log(q/q0)
    bmod_c = b0/r_major*sqrt(r0**2 + r**2/q**2)
    b_theta_c = b0*r**2/(r_major*q)
    b_phi_c = b0*r0
    call prove_circular()
    call emit_circular()


    ! ------------------------------------------------------ torque and weights
    bmod = sym(arena, "bmod")
    call prove_torque()
    call prove_bracket()

    h0 = sym(arena, "h0")
    dens = sym(arena, "dens")
    ddens = sym(arena, "ddens_dpsi")
    temp = sym(arena, "temp")
    dtemp = sym(arena, "dtemp_dpsi")
    phi0 = sym(arena, "phi0")
    dphi0 = sym(arena, "dphi0_dpsi")
    p_phi = sym(arena, "p_phi")
    n_f = func("n_profile", [c_light*p_phi/charge])
    t_f = func("T_profile", [c_light*p_phi/charge])
    phi0_f = func("Phi0_profile", [c_light*p_phi/charge])
    lnf0 = log(n_f) + 3*log(mass/(2*pi_expr(arena)*t_f))/2 &
        - (h0 - charge*phi0_f)/t_f
    dlnf0_closed = c_light/charge*(diff(n_f, p_phi)/n_f*charge/c_light &
        + charge*diff(phi0_f, p_phi)*charge/c_light/t_f &
        + ((h0 - charge*phi0_f)/t_f - num(arena, 3)/2) &
        *diff(t_f, p_phi)*charge/c_light/t_f)
    dlnf0 = diff(lnf0, p_phi)
    ! [P5] chain rule of the canonical Maxwellian against the closed form.
    call check_identity(suite, native, "P5 dlnf0/dp_phi closed form", &
        dlnf0 - dlnf0_closed)
    lnf0 = profile_values(lnf0)
    dlnf0 = profile_values(dlnf0)

    db1_dphi = sym(arena, "db1_dphi")
    dphi1_dphi = sym(arena, "dphi1_dphi")
    dh1_dphi = (mu + mass*vpar**2/bmod)*db1_dphi + charge*dphi1_dphi
    call emit_kernel_file("deltaf_weight_kernel", "evaluate_deltaf_weight", &
        [character(len=16) :: "mass", "charge", "c_light", "mu", "vpar", &
        "bmod", "h0", "dens", "ddens_dpsi", "temp", "dtemp_dpsi", "phi0", &
        "dphi0_dpsi", "db1_dphi", "dphi1_dphi"], &
        [simplified(lnf0), simplified(dlnf0), simplified(dh1_dphi), &
        simplified(dlnf0*dh1_dphi), simplified(-dh1_dphi)], &
        [character(len=16) :: "ln_f0", "dlnf0_dpphi", "dh1_dphi", "dw_dt", &
        "torque_sample"])

    ! ------------------------------------------- collisional background drive
    bmin = sym(arena, "bmin")
    bmax = sym(arena, "bmax")
    branch = sym(arena, "branch")
    dchi = sym(arena, "dchi_dlambda")
    h_phi = sym(arena, "h_phi")
    r_maj = sym(arena, "r_major")
    kappa_f = sym(arena, "dlnf0_dpphi")
    call emit_fax()
    call emit_fax_chartmap()

    call emit_energy_kernels()
    call suite_end(suite)

contains

    !> Harmonic delta B, delta Phi and their (s, theta, phi) derivatives.
    subroutine emit_harmonic()
        type(expr_t) :: roots(8)

        roots(1) = db1
        roots(2) = diff(db1, s)
        roots(3) = diff(db1, theta)
        roots(4) = diff(db1, phi)
        roots(5) = dphi1
        roots(6) = diff(dphi1, s)
        roots(7) = diff(dphi1, theta)
        roots(8) = diff(dphi1, phi)
        roots = harmonic_values(roots)
        call emit_kernel_file("deltaf_harmonic_kernel", &
            "evaluate_deltaf_harmonic", &
            [character(len=16) :: "m_pol", "n_tor", "bmn", "dbmn_ds", "phimn", &
            "dphimn_ds", "theta", "phi"], simplified_all(roots), &
            [character(len=16) :: "db1", "db1_ds", "db1_dtheta", "db1_dphi", &
            "dphi1", "dphi1_ds", "dphi1_dtheta", "dphi1_dphi"])
    end subroutine emit_harmonic

    !> Both real Fourier coefficients, with only the requested derivatives.
    subroutine emit_complex_harmonic()
        type(expr_t) :: bc, bs, field, roots(2), dc, dt

        bc = sym(arena, "bc")
        bs = sym(arena, "bs")
        field = bc*cos(alpha) + bs*sin(alpha)
        dc = n_tor*(bc*sin(alpha) - bs*cos(alpha))
        dt = m_pol*(-bc*sin(alpha) + bs*cos(alpha))
        call check_identity(suite, native, "complex harmonic phi derivative", &
            diff(field, phi) - dc)
        call check_identity(suite, native, "complex harmonic theta derivative", &
            diff(field, theta) - dt)
        roots = [dc, dt]
        call emit_kernel_file("deltaf_complex_harmonic_kernel", &
            "evaluate_deltaf_complex_harmonic", &
            [character(len=16) :: "m_pol", "n_tor", "bc", "bs", "theta", "phi"], &
            simplified_all(roots), &
            [character(len=16) :: "db_dphi", "db_dtheta"])
    end subroutine emit_complex_harmonic

    !> Circular-chart equilibrium.
    subroutine emit_circular()
        type(expr_t) :: roots(11), sv, dpsi_ds

        sv = (r/a_minor)**2
        dpsi_ds = diff(psi_pol, r)/diff(sv, r)
        roots(1) = bmod_c
        roots(2) = diff(bmod_c, r)
        roots(3) = diff(bmod_c, theta)
        roots(4) = psi_pol
        roots(5) = diff(psi_pol, r)
        roots(6) = dpsi_ds
        roots(7) = b_theta_c/bmod_c
        roots(8) = b_phi_c/bmod_c
        roots(9) = sv
        roots(10) = diff(sv, r)
        roots(11) = r_major
        call emit_kernel_file("deltaf_circular_chart_kernel", &
            "evaluate_deltaf_circular_chart", &
            [character(len=16) :: "R0", "a_minor", "B0", "q0", "q1", "r", "theta"], &
            simplified_all(roots), &
            [character(len=16) :: "bmod", "dbmod_dr", "dbmod_dtheta", "psi_pol", &
            "dpsi_dr", "dpsi_ds", "h_theta", "h_phi", "s", "ds_dr", "R_major"])
    end subroutine emit_circular

    !> Coordinates and weight source of the collisional background drive.
    subroutine emit_fax()
        type(expr_t) :: v, lam, mu_lam, gpar, lam_dot, chi_f, hf0, src
        type(expr_t) :: speed, lam_s, src_g
        type(expr_t) :: roots(2)

        ! Symbols of the (x, v, lambda) chart at fixed position.
        v = sym(arena, "v")
        lam = sym(arena, "lambda")
        gpar = sym(arena, "gradpar_b")
        mu_lam = mass*v**2*(1 - lam**2)/(2*bmod)
        ! [F1] lambda' = v_par'/v (v conserved), m v_par' = -mu b.grad B.
        lam_dot = -(mu_lam/(mass*v))*gpar
        call check_identity(suite, native, &
            "F1 dlambda/dt = -(v (1-lambda**2)/2) b.grad ln B", &
            lam_dot + v*(1 - lam**2)/2*gpar/bmod)

        ! Source -delta(lambda') d(h/f0)/dlambda with
        ! delta(lambda') = -(v (1-lambda**2)/2) b.grad(dB/B),
        ! b.grad(dB/B) = (h_phi/R**2) (d dB/dphi)/B for an m = 0 harmonic.
        chi_f = func("chi", [lam])
        hf0 = kappa_f*mass*v*h_phi*chi_f
        src = (v*(1 - lam**2)/2)*(h_phi/r_maj**2)*(db1_dphi/bmod)*diff(hf0, lam)
        ! [F2] the table derivative is the lambda derivative of chi.
        call check_identity(suite, native, "F2 d(h/f0)/dlambda = kappa m v h_phi chi'", &
            diff(hf0, lam) - kappa_f*mass*v*h_phi*partial_derivative(chi_f, 1))
        src = subs_many(src, [partial_derivative(chi_f, 1)], [dchi])
        ! In the caller's variables: v**2 = vpar**2 + 2 mu B/m, lambda = vpar/v.
        speed = sqrt(vpar**2 + 2*mu*bmod/mass)
        lam_s = vpar/speed
        src = subs_many(subs_many(src, [lam], [lam_s]), [v], [speed])

        ! [F3] h = G: chi = lambda, d chi/d lambda = 1; the source is the mu part
        ! of the canonical drive with b^phi h_phi = (h_phi/R)**2.
        src_g = subs_many(src, [dchi], [num(arena, 1)])
        call check_identity(suite, symeng, &
            "F3 chi = lambda gives (dlnf0/dp_phi) mu (h_phi/R)**2 d dB/dphi", &
            src_g - kappa_f*mu*(h_phi/r_maj)**2*db1_dphi)

        roots(1) = lam_s
        roots(2) = branch*acos(1 - 2*(bmod - bmin)/(bmax - bmin))
        call emit_kernel_file("deltaf_fax_coordinates_kernel", &
            "evaluate_deltaf_fax_coordinates", &
            [character(len=16) :: "mass", "mu", "vpar", "bmod", "bmin", "bmax", &
            "branch"], simplified_all(roots), &
            [character(len=16) :: "lambda", "vartheta"])
        call emit_kernel_file("deltaf_fax_weight_kernel", &
            "evaluate_deltaf_fax_weight", &
            [character(len=16) :: "mu", "h_phi", "r_major", "dlnf0_dpphi", &
            "db1_dphi", "dchi_dlambda"], &
            [simplified(src)], [character(len=16) :: "dw_dt_fax"])
    end subroutine emit_fax

    !> Fixed covariant Boozer B_k: h^i scales with |B|, so variation of
    !> h^i d_i ln B is h^i d_i(delta B)/B (the denominator terms cancel).
    !> Thus -delta(lambda_dot) d(h/f0)/dlambda = kappa mu h_phi h^i d_i B1 chi'.
    subroutine emit_fax_chartmap()
        type(expr_t) :: ht, ct, cp, gt, gp, dbt, rb1, denom, src, gradpar

        ht = sym(arena, "h_theta")
        ct = sym(arena, "curl_theta")
        cp = sym(arena, "curl_phi")
        gt = sym(arena, "gb_theta")
        gp = sym(arena, "gb_phi")
        dbt = sym(arena, "db1_dtheta")
        rb1 = sym(arena, "relative_db1")
        denom = ht*ct + h_phi*cp
        ! rb1=0 for absolute amplitudes; otherwise include derivatives of
        ! the local B multiplying the relative harmonic amplitude.
        src = kappa_f*mu*h_phi*(ct*(dbt + rb1*gt/bmod) &
            + cp*(db1_dphi + rb1*gp/bmod))/denom*dchi
        gradpar = (ct*gt + cp*gp)/denom
        call emit_kernel_file("deltaf_fax_chartmap_kernel", &
            "evaluate_deltaf_fax_chartmap", &
            [character(len=16) :: "mu", "h_theta", "h_phi", "curl_theta", &
            "curl_phi", "gb_theta", "gb_phi", "dlnf0_dpphi", &
            "db1_dtheta", "db1_dphi", "dchi_dlambda", "bmod", "relative_db1"], &
            simplified_all([src, gradpar]), &
            [character(len=16) :: "dw_dt_fax", "gradpar_b"])
    end subroutine emit_fax_chartmap

    !> Circular-chart identities.
    subroutine prove_circular()
        ! GORILLA field_analytic_circ_mod: B_Z = (1/R) dpsi/dR = B0 (R-R0)/(R q)
        ! along r at theta = 0, i.e. dpsi/dr = B0 r/q.
        call check_identity(suite, native, "C1 dpsi/dr = B0 r/q", &
            diff(psi_pol, r) - b0*r/q)
        ! Unit field: B**2 = (B_phi/R)**2 + (B_theta/r)**2 for the covariant
        ! components of the orthogonal (r, theta, phi) chart.
        call check_identity(suite, symeng, "C2 |h| = 1", &
            (b_phi_c/r_major)**2 + (b_theta_c/r)**2 - bmod_c**2)
    end subroutine prove_circular

    !> [P1], [P2]: the canonical torque from the Euler-Lagrange equation.
    subroutine prove_torque()
        type(expr_t) :: b_s, b_t, b_p, bm, xd_s, xd_t, xd_p, el_force, kappa
        type(expr_t) :: h_c, phi_f, sv

        sv = sym(arena, "s")
        bm = func("B", [sv, theta, phi])
        b_s = func("B_s", [sv, theta])
        b_t = func("B_theta", [sv, theta])
        b_p = func("B_phi", [sv, theta])
        xd_s = sym(arena, "sdot")
        xd_t = sym(arena, "thetadot")
        ! The v_par Euler-Lagrange equation h_k x'^k = v_par fixes phidot.
        xd_p = (vpar*bm - b_s*xd_s - b_t*xd_t)/b_p
        el_force = mass*vpar*(diff(b_s/bm, phi)*xd_s + diff(b_t/bm, phi)*xd_t &
            + diff(b_p/bm, phi)*xd_p)
        call check_identity(suite, native, &
            "P1 m v_par dh_k/dphi x'^k = -m v_par**2/B dB/dphi", &
            el_force + mass*vpar**2/bm*diff(bm, phi))

        ! Canonical H at fixed (s, kappa = v_par/B), determined by p_theta and
        ! p_phi without reference to phi for flux-function B_theta, B_phi.
        kappa = sym(arena, "kappa")
        phi_f = func("Phi", [sv, theta, phi])
        h_c = mass*kappa**2*bm**2/2 + mu*bm + charge*phi_f
        call check_identity(suite, native, &
            "P2 dH/dphi|canonical = (mu + m v_par**2/B) dB/dphi + e dPhi/dphi", &
            subs_many(diff(h_c, phi), [kappa], [vpar/bm]) &
            - ((mu + mass*vpar**2/bm)*diff(bm, phi) + charge*diff(phi_f, phi)))
    end subroutine prove_torque

    !> [P3], [P4]: stationarity of f0 and the linearised Vlasov source.
    subroutine prove_bracket()
        type(expr_t) :: hz, h1, f0

        th = sym(arena, "q_theta")
        pth = sym(arena, "p_theta")
        ph = sym(arena, "q_phi")
        pph = sym(arena, "p_phi_c")
        hz = func("H0", [th, pth, pph])
        h1 = func("H1", [th, pth, ph, pph])
        f0 = func("F", [pph, hz])
        call check_identity(suite, native, "P3 {f0, H0} = 0", &
            bracket(f0, hz))
        ! -{f0,H1} = dF/dp_phi dH1/dphi + dF/dH0 {H1, H0}.
        call check_identity(suite, native, &
            "P4 -{f0,H1} = F_p dH1/dphi + F_H dH1/dt|0", &
            -bracket(f0, h1) - (fp(pph, hz)*diff(h1, ph) &
            + fh(pph, hz)*bracket(h1, hz)))
    end subroutine prove_bracket

    !> Canonical bracket in (q_theta, p_theta), (q_phi, p_phi_c).
    function bracket(fa, fb) result(res)
        type(expr_t), intent(in) :: fa, fb
        type(expr_t) :: res

        res = diff(fa, th)*diff(fb, pth) - diff(fa, pth)*diff(fb, th) &
            + diff(fa, ph)*diff(fb, pph) - diff(fa, pph)*diff(fb, ph)
    end function bracket

    !> dF/dp_phi at fixed H0 and dF/dH0, as FortSym's canonical partials.
    function fp(pv, hv) result(res)
        type(expr_t), intent(in) :: pv, hv
        type(expr_t) :: res
        type(expr_t) :: x

        x = sym(arena, "fp_probe")
        res = subs_many(diff(func("F", [x, hv]), x), [x], [pv])
    end function fp

    function fh(pv, hv) result(res)
        type(expr_t), intent(in) :: pv, hv
        type(expr_t) :: res
        type(expr_t) :: x

        x = sym(arena, "fh_probe")
        res = subs_many(diff(func("F", [pv, x]), x), [x], [hv])
    end function fh

    !> Replace the opaque harmonic profiles by caller data symbols.
    function harmonic_values(e) result(res)
        type(expr_t), intent(in) :: e(:)
        type(expr_t) :: res(size(e))
        integer :: k

        do k = 1, size(e)
            res(k) = subs_many(e(k), [diff(b_of_s, s), diff(phi_of_s, s), &
                b_of_s, phi_of_s], [dbmn_ds, dphimn_ds, bmn, phimn])
        end do
    end function harmonic_values

    !> Replace the opaque profiles n, T, Phi0 at psi* by caller data symbols.
    function profile_values(e) result(res)
        type(expr_t), intent(in) :: e
        type(expr_t) :: res

        res = subs_many(e, [partial_derivative(n_f, 1), &
            partial_derivative(t_f, 1), partial_derivative(phi0_f, 1), &
            n_f, t_f, phi0_f], [ddens, dtemp, dphi0, dens, temp, phi0])
    end function profile_values

    ! Horner step and physical f_AX/f_M -> remainder h/f_M. I/O and loops
    ! remain handwritten; these are the new pointwise mathematical leaves.
    subroutine emit_energy_kernels()
        type(expr_t) :: xx, aa, cc, gg, vv, hh, kk, mm, bb, ee, tt
        xx = sym(arena, "x")
        aa = sym(arena, "acc")
        cc = sym(arena, "coefficient")
        call emit_kernel_file("deltaf_fax_horner_kernel", "evaluate_fax_horner", &
            [character(len=16) :: "x", "acc", "coefficient"], &
            [xx*aa+cc], [character(len=16) :: "value"])
        gg = sym(arena, "dg")
        vv = sym(arena, "speed")
        hh = sym(arena, "hphi")
        kk = sym(arena, "kappa")
        mm = sym(arena, "mass")
        ee = sym(arena, "energy")
        tt = sym(arena, "temperature")
        bb = gg/(mm*vv*hh)-kk
        call emit_kernel_file("deltaf_fax_energy_kernel", "evaluate_fax_energy", &
            [character(len=16) :: "dg", "mass", "speed", "hphi", "kappa"], &
            [bb], [character(len=16) :: "effective_kappa"])
        call emit_kernel_file("deltaf_fax_speed_kernel", "evaluate_fax_speed", &
            [character(len=16) :: "energy", "mass", "temperature"], &
            [sqrt(2*ee/mm), sqrt(ee/tt)], &
            [character(len=16) :: "speed", "x"])
        call check_identity(suite, symeng, "f_AX=G removes remainder", &
            simplified(bb*mm*vv*hh - (gg-kk*mm*vv*hh)))
    end subroutine emit_energy_kernels

    function simplified(e) result(res)
        type(expr_t), intent(in) :: e
        type(expr_t) :: res
        type(engine_result_t) :: rr

        rr = symeng%simplify(e)
        if (rr%ok) then
            res = rr%value
        else
            res = e
        end if
    end function simplified

    function simplified_all(e) result(res)
        type(expr_t), intent(in) :: e(:)
        type(expr_t) :: res(size(e))
        integer :: k

        do k = 1, size(e)
            res(k) = simplified(e(k))
        end do
    end function simplified_all

    subroutine emit_kernel_file(module_name, procedure_name, arg_names, roots, &
            output_names)
        character(*), intent(in) :: module_name, procedure_name
        character(*), intent(in) :: arg_names(:), output_names(:)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t) :: spec
        type(str_t) :: emitted
        character(:), allocatable :: path
        integer :: unit, ios, k
        logical :: ok

        if (size(roots) /= size(output_names)) error stop "arity mismatch"
        spec%name = str(procedure_name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%temp_prefix = str("t")
        spec%generator = str(GENERATOR)
        spec%generator_revision = str(FORTSYM_REVISION)
        spec%regenerate_command = str(REGENERATE)
        allocate (spec%args(size(arg_names)), spec%outputs(size(output_names)))
        do k = 1, size(arg_names)
            spec%args(k) = str(trim(arg_names(k)))
        end do
        do k = 1, size(output_names)
            spec%outputs(k) = str(trim(output_names(k)))
        end do
        emitted = emit_kernel(roots, spec, ok)
        if (.not. ok .or. len(chars(emitted)) == 0) then
            write (output_unit, "(a)") "refused: "//module_name
            do k = 1, size(roots)
                write (output_unit, "(a)") "  root: "//chars(print_expr(roots(k)))
            end do
            error stop "fortsym refused a generated kernel"
        end if
        path = trim(output_dir)//"/"//module_name//".f90"
        open (newunit=unit, file=path, status="replace", action="write", &
            iostat=ios)
        if (ios /= 0) error stop "cannot open generated kernel output"
        write (unit, "(a)") chars(emitted)
        close (unit)
        write (output_unit, "(a)") "wrote "//path
    end subroutine emit_kernel_file

end program gen_deltaf_kernels
