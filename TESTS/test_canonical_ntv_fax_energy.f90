program test_canonical_ntv_fax_energy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use deltaf_fax_energy_kernel, only: evaluate_fax_energy
    use canonical_ntv_mod, only: ntv_fax_load_table, ntv_fax_response
    implicit none
    real(dp) :: mass(2), charge(2), temperature(2), amplitude(3,2)
    real(dp) :: coeff(0:2,2,1), edges(2), centers(1), vt(3), lam(3)
    real(dp) :: values(3,3,2,3,2,2), derivatives(3,3,2,3,2,2)
    real(dp) :: energy, x, theta, pitch, expected, got, err
    integer :: u, i, j, m, k, b, a, t, ios
    character(len=512) :: path, arg

    if (command_argument_count() == 2) then
        call get_command_argument(1, path)
        call ntv_fax_load_table(trim(path))
        call get_command_argument(2, path)
        open(newunit=u,file=trim(path),status='old')
        do
            read(u,*,iostat=ios) a,theta,pitch,energy
            if (ios < 0) exit
            if (ios /= 0) error stop 'invalid oracle sample'
            print '(a, 2es26.17)', 'RESPONSE ', &
                ntv_fax_response(theta,pitch,energy,a,.false.), &
                ntv_fax_response(theta,pitch,energy,a,.true.)
        end do
        close(u)
        stop
    end if

    ! Also serves as the native consumer for the cross-code sample oracle.
    if (command_argument_count() == 5) then
        call get_command_argument(1, path)
        call get_command_argument(2, arg); read(arg,*) a
        call get_command_argument(3, arg); read(arg,*) theta
        call get_command_argument(4, arg); read(arg,*) pitch
        call get_command_argument(5, arg); read(arg,*) energy
        call ntv_fax_load_table(trim(path))
        print '(a, 2es26.17)', 'RESPONSE ', &
            ntv_fax_response(theta,pitch,energy,a,.false.), &
            ntv_fax_response(theta,pitch,energy,a,.true.)
        stop
    end if
    ! Numerical cancellation and a zero profile gradient must both be regular.
    call evaluate_fax_energy(0.125_dp,2.0_dp,0.5_dp,0.25_dp,0.5_dp,got)
    if (abs(got) > 1.0e-15_dp) error stop 'f_AX=G canonical shift'
    call evaluate_fax_energy(0.125_dp,2.0_dp,0.5_dp,0.25_dp,0.0_dp,got)
    if (abs(got-0.5_dp) > 1.0e-15_dp) error stop 'zero profile gradient'
    mass = [1.0_dp, 3.0_dp]
    charge = [-1.0_dp, 1.0_dp]
    temperature = [2.0_dp, 5.0_dp]
    amplitude(:,1) = [0.3_dp, -0.4_dp, 0.8_dp]
    amplitude(:,2) = [0.9_dp, 0.5_dp, -0.2_dp]
    coeff = 0.0_dp
    coeff(0,1,1) = 1.0_dp
    coeff(:,2,1) = [2.5_dp, 0.0_dp, -1.0_dp]
    edges = [0.0_dp, 4.0_dp]
    centers = 0.0_dp
    vt = [-acos(-1.0_dp), 0.0_dp, acos(-1.0_dp)]
    lam = [-1.0_dp, 0.0_dp, 1.0_dp]
    do a = 1, 2
        do b = 1, 2
            do k = 1, 3
                do m = 1, 2
                    do j = 1, 3
                        do i = 1, 3
                            values(i,j,m,k,b,a) = a+2*b+3*k+4*m &
                                + vt(j)*lam(i)
                            derivatives(i,j,m,k,b,a) = vt(j)
                        end do
                    end do
                end do
            end do
        end do
    end do
    open(newunit=u,file='fax_energy_test.dat',status='replace')
    write(u,'(a)') 'FAX_ENERGY_V1'
    write(u,*) 3,3,2,2,2,1,.true.
    write(u,*) mass,charge,temperature,amplitude,edges,centers,coeff,vt,lam, &
        values,derivatives
    close(u)
    call ntv_fax_load_table('fax_energy_test.dat')
    err = 0.0_dp
    do a = 1, 2
        do t = 0, 16
            x = 0.3_dp*t
            energy = temperature(a)*x*x
            theta = 0.1_dp*t
            pitch = -0.8_dp+0.1_dp*t
            expected = 0.0_dp
            do b = 1, 2
                do k = 1, 3
                    expected = expected+amplitude(k,b) * &
                        ((a+2*b+3*k+4+theta*pitch) &
                        +(2.5_dp-x*x)*(a+2*b+3*k+8+theta*pitch))
                end do
            end do
            got = ntv_fax_response(theta,pitch,energy,a,.false.)
            err = max(err,abs(got-expected)/max(1.0_dp,abs(expected)))
            expected = sum(amplitude)*theta*(3.5_dp-x*x)
            got = ntv_fax_response(theta,pitch,energy,a,.true.)
            err = max(err,abs(got-expected)/max(1.0_dp,abs(expected)))
        end do
    end do
    print '(a, es24.16)', 'energy/species max normalized error: ', err
    if (err > 1.0e-13_dp) error stop 'f_AX energy/species reconstruction'
    print '(a)', 'PASS: f_AX energy/species reconstruction'
end program test_canonical_ntv_fax_energy
