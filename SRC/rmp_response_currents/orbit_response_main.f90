program gorilla_orbit_response
    ! Backward fractional-delta-f characteristics, paired with an exact cylinder
    ! control variate. The finite-R correction uses the actual tetrahedral pusher.
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use omp_lib, only: omp_set_dynamic, omp_get_num_threads, omp_get_wtime
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use cylindrical_electron_response_mod, only: assemble_cylindrical_electrons, close_energy_response,assemble_diffused_cylinder
    use orbit_timestep_gorilla_mod, only: initialize_gorilla, orbit_timestep_gorilla
    use gorilla_settings_mod, only: load_gorilla_inp
    use tetra_grid_settings_mod, only: load_tetra_grid_inp, R0_analytic_circ, n_field_periods_manual, &
        sfc_s_min, sfc_s_max, a_analytic_circ,n3
    use collis_ions, only: exact_ou_velocity
    use supporting_functions_mod, only: bmod_func, energy_tot_func
    use tetra_physics_mod, only: tetra_physics
    use field_analytic_circ_mod, only: psi_pol_analytic_circ
    use find_tetra_mod, only: find_tetra
    use utils_rmp_response_currents_mod, only: draw_annulus_rphiz_analytic
    use anomalous_transport_displacement_mod, only: anomalous_transport_displacement,compute_diffusion_cholesky
    implicit none
    real(dp), parameter :: pi=acos(-1.0_dp)
    complex(dp), parameter :: ii=(0.0_dp,1.0_dp)
    character(1024) :: input_path,output_path,header,arg
    integer :: unit,M,N,model,mmode,nmode,j,k,col,d,ir,im,seed_size,ib,phase_unit,channels
    integer :: radial_samples=64,markers=32,batches=4,seed=18371,refinement=1
    integer :: n_respawn_max=0,n_respawn_total=0,n_respawn_markers=0
    integer :: n_respawn_inside=0,n_respawn_outside=0,max_respawns_per_marker=0,respawn_unit
    real(dp) :: lag_nu=8.0_dp,nu_dt=0.05_dp,phase_step=0.15_dp
    integer :: orbit_threads=1
    logical :: cylinder_only=.false.,marker_streams=.false.
    real(dp) :: anomalous_diffusion_coefficient=0.0_dp, diffusion_step_cm=0.1_dp
    integer :: transport_boundary=0
    integer :: diffusion_kicks=0,diffusion_reflections=0,diffusion_respawns=0
    real(dp) :: diffusion_bounds(2)
    real(dp) :: diffusion_edge_flux
    namelist /orbit_response/ radial_samples,markers,batches,seed,lag_nu,nu_dt,phase_step,cylinder_only,n_respawn_max, &
        orbit_threads,marker_streams,anomalous_diffusion_coefficient,diffusion_step_cm,transport_boundary
    real(dp) :: L,rm,clight,r0,bg0(13),qn,coef0,coef,end_t
    real(dp) :: closure_error,rcond,max_closure_error=0.0_dp,min_closure_rcond=1.0_dp
    real(dp) :: max_radial_excursion=0.0_dp,phase_error=0.0_dp,spawn_error=0.0_dp
    real(dp) :: max_energy_error=0.0_dp,max_reversal_error=0.0_dp,zero_v_phase_error=0.0_dp
    integer,parameter :: flux_points=8193
    real(dp) :: flux_r(flux_points),flux_psi(flux_points)
    real(dp),allocatable :: background(:,:),wave(:),seed_draws(:,:)
    integer,allocatable :: seed_array(:),marker_seeds(:,:),master_seed(:)
    complex(dp),allocatable :: blocks(:,:,:),corrections(:,:,:,:),marker_local(:,:,:)
    complex(dp),allocatable :: batch_local(:,:,:),output_phase(:)
    complex(dp),allocatable :: bare_blocks(:,:,:),cylinder_closed(:,:,:),closed(:,:,:),mean_blocks(:,:,:)
    complex(dp),allocatable :: diffused_reference(:,:,:)
    type marker_report
        real(dp) :: max_radial_excursion=0.0_dp,phase_error=0.0_dp,spawn_error=0.0_dp
        real(dp) :: max_energy_error=0.0_dp,max_reversal_error=0.0_dp,zero_v_phase_error=0.0_dp
        real(dp) :: phase_row(5)=0.0_dp
        integer :: n_respawn_total=0,n_respawn_markers=0,n_respawn_inside=0,n_respawn_outside=0
        integer :: max_respawns_per_marker=0
        integer :: diffusion_kicks=0,diffusion_reflections=0,diffusion_respawns=0
        real(dp),allocatable :: events(:,:)
    end type
    type,extends(marker_report) :: marker_state
        real(dp) :: theta0,alpha0,alpha,dt,t,vp
        real(dp) :: vperp,vc,v0,vc_old,vc_phase,p0
        real(dp) :: zeta,rnow,coef,phase_rate,predict_phase,phase_probe
        real(dp) :: probe_dt,geometric_radius,probe_vp,probe_vperp,probe_mu,energy0
        real(dp) :: energy1
        real(dp) :: reference_radial_displacement=0,reference_diffusion_phase=0
        real(dp) :: reference_factor(3,3),reference_radial_covector(3),reference_phase_covector(3)
        real(dp) :: x(3)
        real(dp) :: bgt(13)
        real(dp) :: rand2(2)
        real(dp) :: probe_x(3)
        real(dp) :: delta_x(3)
        integer :: tetr,iface,probe_tetr,probe_iface
        integer :: n_respawn_used
        logical :: initialized,probe_initialized
        complex(dp),allocatable :: local(:,:),previous(:,:),present_source(:,:),source_phase(:)
        complex(dp) :: hphase,cphase
    end type
    type(marker_report),allocatable :: reports(:)
    integer :: event,actual_threads=1,progress_unit,timing_unit
    real(dp) :: mesh_start,mesh_seconds,trace_start,trace_seconds,closure_start
    call get_command_argument(1,input_path)
    call get_command_argument(2,output_path)
    call get_command_argument(3,arg)
    if(len_trim(arg)>0) read(arg,*) refinement
    open(newunit=unit,file=trim(input_path),status='old')
    read(unit,'(a)') header
    if(trim(header)/='GK_BACKGROUND_V1') error stop 'wrong background format'
    read(unit,*) M,N,model,mmode,nmode
    read(unit,*) L,rm,clight
    if(model/=0.and.model/=1) error stop 'orbit provider supports electron models 0/1'
    allocate(background(N,13))
    do j=1,N
        read(unit,*) background(j,:)
    end do
    close(unit)
    open(newunit=unit,file='orbit_response.nml',status='old')
    read(unit,nml=orbit_response)
    close(unit)
    if(radial_samples<2*M+1.or.markers<1.or.batches<2.or.mod(markers,batches)/=0) &
        error stop 'need radial_samples >= 2M+1 and equal independent batches'
    if(min(lag_nu,nu_dt,phase_step)<=0.0_dp) error stop 'invalid orbit quadrature'
    if(orbit_threads<1) error stop 'orbit_threads must be positive'
    if(orbit_threads>1.and..not.marker_streams) error stop 'threading requires marker_streams'
    call omp_set_dynamic(.false.)
    if(n_respawn_max<0) error stop 'negative respawn budget'
    if(.not.ieee_is_finite(anomalous_diffusion_coefficient).or.anomalous_diffusion_coefficient<0.0_dp) &
        error stop 'invalid anomalous diffusion coefficient'
    if(anomalous_diffusion_coefficient>0.0_dp) then
        if(cylinder_only) error stop 'spatial diffusion requires actual mesh trajectories'
        if(transport_boundary/=1.and.transport_boundary/=2) error stop 'positive diffusion requires explicit radial boundary'
        if(transport_boundary==2.and.n_respawn_max==0) error stop 'diffusion respawn requires a positive respawn budget'
        if(.not.ieee_is_finite(diffusion_step_cm).or.diffusion_step_cm<=0.0_dp) &
            error stop 'invalid diffusion substep length'
    end if
    call random_seed(size=seed_size)
    allocate(seed_array(seed_size))
    seed_array=seed
    call random_seed(put=seed_array)
    call load_gorilla_inp()
    call load_tetra_grid_inp()
    mesh_start=omp_get_wtime()
    call initialize_gorilla(ipert_in=0)
    if(anomalous_diffusion_coefficient>0.0_dp.and.transport_boundary==1) then
        ! grid_kind=5 maps s to geometric rho through normalized toroidal flux.
        ! Its outer boundary is a regular n3-gon in (R-R0,Z), not a circle.
        ! The apothem is the largest circular reflecting surface inside it.
        diffusion_edge_flux=a_analytic_circ**2/(R0_analytic_circ+sqrt(R0_analytic_circ**2-a_analytic_circ**2))
        diffusion_bounds=sqrt([sfc_s_min,sfc_s_max]*diffusion_edge_flux &
            *(2*R0_analytic_circ-[sfc_s_min,sfc_s_max]*diffusion_edge_flux))
        diffusion_bounds(2)=diffusion_bounds(2)*cos(pi/n3)
        if(diffusion_bounds(1)>=diffusion_bounds(2)) error stop 'mesh-contained reflecting annulus is empty'
    end if
    mesh_seconds=omp_get_wtime()-mesh_start
    do j=1,flux_points
        flux_r(j)=1.0_dp+49.0_dp*real(j-1,dp)/real(flux_points-1,dp)
        flux_psi(j)=abs(psi_pol_analytic_circ(flux_r(j)))
    end do
    d=2*M+1
    channels=4
    if(model==1) then
        channels=9
        call assemble_cylindrical_electrons(background,M,L,clight,blocks,refinement,model,bare_blocks)
        call close_energy_response(bare_blocks,cylinder_closed,closure_error,rcond)
        max_closure_error=max(max_closure_error,closure_error)
        min_closure_rcond=min(min_closure_rcond,rcond)
    else
        call assemble_cylindrical_electrons(background,M,L,clight,blocks,refinement)
    end if
    if(anomalous_diffusion_coefficient>0) then
        ! Preserve the laminar model1 projection correction exactly. Replace
        ! only the BARE traced control variate with its diffused expectation.
        if(model==1) call write_response('orbit_energy_base.dat',bare_blocks)
        call assemble_diffused_cylinder(background,M,L,clight,anomalous_diffusion_coefficient,channels,diffused_reference)
        call write_response('orbit_transport_reference.dat',diffused_reference)
        if(model==1) then
            call move_alloc(diffused_reference,bare_blocks)
        else
            call move_alloc(diffused_reference,blocks)
        end if
    end if
    allocate(wave(d),corrections(d,d,channels,batches),marker_local(d,9,markers),reports(markers))
    if(marker_streams) allocate(marker_seeds(seed_size,markers),seed_draws(seed_size,markers),master_seed(seed_size))
    allocate(batch_local(d,channels,batches))
    allocate(output_phase(d))
    wave=2*pi*real([(j-M-1,j=1,d)],dp)/L
    corrections=0
    open(newunit=phase_unit,file="orbit_phase_samples.dat",status="new")
    write(phase_unit,'(a)') "# r_cm u0 theta0 actual_frequency_over_nu KIM_frequency_over_nu"
    if(n_respawn_max>0) then
        open(newunit=respawn_unit,file='orbit_respawns.dat',status='new')
        write(respawn_unit,'(a)') '# sample marker lag_s pending_step_s old_geometric_r_cm new_geometric_r_cm ' &
            //'old_vpar old_vperp new_vpar new_vperp'
    end if
    open(newunit=progress_unit,file='orbit_progress.dat',status='new')
    write(progress_unit,'(a)') '# radial_sample cumulative_trace_seconds'
    trace_start=omp_get_wtime()
    do ir=1,radial_samples
        ! Midpoint quadrature of the same absolute, periodized radial coordinate.
        r0=background(1,1)+(real(ir,dp)-0.5_dp)*L/radial_samples
        call sample_background(r0,bg0)
        qn=bg0(5)**(-2)/(4*pi)*bg0(7)**2*bg0(13)/(bg0(8)*clight)
        output_phase=exp(-ii*wave*r0)/real(radial_samples,dp)
        coef0=bg0(9)+bg0(10)
        end_t=lag_nu/bg0(6)
        batch_local=0
        if(marker_streams) then
            ! Reserve all stream seeds serially, independent of adaptive path length.
            ! Restore the reservation stream after workers seed/consume their own RNGs.
            call random_number(seed_draws)
            marker_seeds=int(seed_draws*real(huge(1),dp))
            call random_seed(get=master_seed)
            !$omp parallel default(none) num_threads(orbit_threads) &
            !$omp shared(marker_seeds,marker_local,reports,ir,markers,actual_threads)
            !$omp single
            actual_threads=omp_get_num_threads()
            !$omp end single
            !$omp do schedule(dynamic,1)
            do im=1,markers
                call trace_marker(ir,im,marker_local(:,:,im),reports(im),marker_seeds(:,im))
            end do
            !$omp end do
            !$omp end parallel
            call random_seed(put=master_seed)
        else
            ! Default path retains the original single serial random realization.
            do im=1,markers
                call trace_marker(ir,im,marker_local(:,:,im),reports(im))
            end do
        end if
        ! Deterministic reduction and log order: no worker writes shared files
        ! or adds to shared matrices, including on soft-respawn paths.
        do im=1,markers
            ib=mod(im-1,batches)+1
            batch_local(:,:,ib)=batch_local(:,:,ib)+marker_local(:,1:channels,im)
            call reduce_report(reports(im))
            if(im==1) write(phase_unit,'(5es24.16)') reports(im)%phase_row
            do event=1,reports(im)%n_respawn_total
                write(respawn_unit,'(2i10,8es24.16)') ir,im,reports(im)%events(:,event)
            end do
        end do
        if(n_respawn_max>0) flush(respawn_unit)
        ! At fixed output radius, output_phase and q*n are marker independent.
        ! Sum bare source vectors within each original Monte Carlo batch first;
        ! form one outer product per batch, before the nonlinear energy closure.
        ! Deposition retains the batch estimator for either explicit RNG policy.
        do ib=1,batches
            do k=1,channels
                ! Energy feedback channels 5,6,9 are moments of fractional g.
                coef=qn
                if(k==5.or.k==6.or.k==9) coef=1.0_dp
                do col=1,d
                    corrections(:,col,k,ib)=corrections(:,col,k,ib) &
                        +output_phase*batch_local(col,k,ib)*coef/real(markers/batches,dp)
                end do
            end do
        end do
        write(progress_unit,'(i8,es24.16)') ir,omp_get_wtime()-trace_start
        flush(progress_unit)
        print *, 'Orbit radial sample ',ir,'/',radial_samples
        flush(6)
    end do
    trace_seconds=omp_get_wtime()-trace_start
    closure_start=omp_get_wtime()
    close(progress_unit)
    close(phase_unit)
    if(n_respawn_max>0) close(respawn_unit)
    if(model==1) then
        ! Close the mean bare operator first. The control variate subtracts
        ! the identically Fourier-projected cylinder closure and restores its
        ! exact local model1 result, never a fitted/native KIM correction.
        call close_energy_response(bare_blocks+sum(corrections,dim=4)/batches,closed,closure_error,rcond)
        max_closure_error=max(max_closure_error,closure_error)
        min_closure_rcond=min(min_closure_rcond,rcond)
        mean_blocks=blocks+closed-cylinder_closed
        if(anomalous_diffusion_coefficient==0) call write_response('orbit_energy_base.dat',bare_blocks)
        call write_response('orbit_energy_exact.dat',blocks)
        call write_response('orbit_energy_mean.dat',bare_blocks+sum(corrections,dim=4)/batches)
    else
        mean_blocks=blocks+sum(corrections,dim=4)/batches
    end if
    call write_response(trim(output_path),mean_blocks)
    do ib=1,batches
        write(arg,'(a,i0,a)') 'orbit_batch_',ib,'.dat'
        if(model==1) then
            deallocate(closed)
            call close_energy_response(bare_blocks+corrections(:,:,:,ib),closed,closure_error,rcond)
            max_closure_error=max(max_closure_error,closure_error)
            min_closure_rcond=min(min_closure_rcond,rcond)
            call write_response(trim(arg),blocks+closed-cylinder_closed)
            write(arg,'(a,i0,a)') 'orbit_energy_batch_',ib,'.dat'
            call write_response(trim(arg),bare_blocks+corrections(:,:,:,ib))
        else
            call write_response(trim(arg),blocks+corrections(:,:,:,ib))
        end if
    end do
    open(newunit=unit,file='orbit_diagnostics.dat',status='new')
    write(unit,'(a,i0)') 'orbit_threads ',actual_threads
    write(unit,'(a,i0)') 'marker_streams ',merge(1,0,marker_streams)
    if(anomalous_diffusion_coefficient>0.0_dp) then
        write(unit,'(a,es24.16)') 'anomalous_diffusion_cm2_per_s ',anomalous_diffusion_coefficient
        write(unit,'(a,i0)') 'transport_boundary ',transport_boundary
        write(unit,'(a,i0)') 'diffusion_respawns ',diffusion_respawns
        write(unit,'(a,es24.16)') 'diffusion_step_cm ',diffusion_step_cm
        write(unit,'(a,i0)') 'diffusion_kicks ',diffusion_kicks
        write(unit,'(a,i0)') 'diffusion_reflections ',diffusion_reflections
        if(transport_boundary==1) then
            write(unit,'(a,es24.16)') 'diffusion_lower_boundary_cm ',diffusion_bounds(1)
            write(unit,'(a,es24.16)') 'diffusion_upper_boundary_cm ',diffusion_bounds(2)
        end if
    end if
    write(unit,'(a,i0)') 'response_markers ',radial_samples*markers
    write(unit,'(a,i0)') 'respawn_events ',n_respawn_total
    write(unit,'(a,i0)') 'respawned_markers ',n_respawn_markers
    write(unit,'(a,i0)') 'respawns_inside_annulus ',n_respawn_inside
    write(unit,'(a,i0)') 'respawns_outside_annulus ',n_respawn_outside
    write(unit,'(a,i0)') 'max_respawns_per_marker ',max_respawns_per_marker
    write(unit,'(a,es24.16)') 'max_collisionless_relative_energy_error ',max_energy_error
    write(unit,'(a,es24.16)') 'max_collisionless_reversal_error ',max_reversal_error
    write(unit,'(a,es24.16)') 'max_flux_radius_spawn_error_cm ',spawn_error
    write(unit,'(a,es24.16)') 'max_radial_excursion_cm ',max_radial_excursion
    write(unit,'(a,es24.16)') 'max_phase_frequency_error_over_nu ',phase_error
    write(unit,'(a,es24.16)') 'max_zero_velocity_phase_frequency_error_over_nu ',zero_v_phase_error
    if(model==1) then
        write(unit,'(a,es24.16)') 'max_energy_closure_backward_residual ',max_closure_error
        write(unit,'(a,es24.16)') 'min_energy_closure_reciprocal_condition_number ',min_closure_rcond
        write(unit,'(a,es24.16)') 'relative_projected_cylinder_closure_difference ', &
            sqrt(sum(abs(cylinder_closed-blocks)**2)/sum(abs(blocks)**2))
    end if
    do k=1,4
        write(unit,'(a,i0,1x,es24.16)') 'relative_block_correction_',k, &
            sqrt(sum(abs(mean_blocks(:,:,k)-blocks(:,:,k))**2)/sum(abs(blocks(:,:,k))**2))
    end do
    close(unit)
    open(newunit=timing_unit,file='orbit_timing.dat',status='new')
    write(timing_unit,'(a,es24.16)') 'mesh_initialization_seconds ',mesh_seconds
    write(timing_unit,'(a,es24.16)') 'tracing_seconds ',trace_seconds
    write(timing_unit,'(a,es24.16)') 'closure_output_seconds ',omp_get_wtime()-closure_start
    close(timing_unit)
contains
    subroutine trace_marker(ir,im,result,report,stream_seed)
        integer,intent(in) :: ir,im
        integer,intent(in),optional :: stream_seed(:)
        complex(dp),intent(out) :: result(d,9)
        type(marker_report),intent(out) :: report
        type(marker_state) :: s
        integer :: spawn_iteration
        allocate(s%local(d,9),s%previous(d,3),s%present_source(d,3),s%source_phase(d))
        if(present(stream_seed)) call random_seed(put=stream_seed)
        s%n_respawn_used=0
        call random_number(s%rand2)
        s%v0=bg0(7)*sqrt(-2*log(max(s%rand2(1),tiny(1.0_dp))))*cos(2*pi*s%rand2(2))
        ! Average perpendicular source energy, retaining thermal mu in the orbit.
        call random_number(s%rand2)
        s%theta0=2*pi*s%rand2(1)
        s%geometric_radius=r0
        s%vp=s%v0
        s%vc=s%v0
        s%vperp=sqrt(2.0_dp)*bg0(7)
        do spawn_iteration=1,5
            s%x=[R0_analytic_circ+s%geometric_radius*cos(s%theta0), &
                2*pi*s%rand2(2)/n_field_periods_manual,s%geometric_radius*sin(s%theta0)]
            s%initialized=.false.
            s%tetr=0
            s%iface=0
            call orbit_timestep_gorilla(s%x,s%vp,s%vperp,0.0_dp,s%initialized,s%tetr,s%iface)
            if(s%tetr==-1) error stop 'initial point outside orbit mesh'
            s%rnow=flux_radius(s%x,s%tetr)
            if(abs(s%rnow-r0)<1e-9_dp) exit
            s%geometric_radius=s%geometric_radius+r0-s%rnow
        end do
        s%spawn_error=max(s%spawn_error,abs(s%rnow-r0))
        if(abs(s%rnow-r0)>1e-8_dp) error stop 'could not match output flux radius'
        s%alpha0=mmode*theta_sfl(s%x)+nmode*n_field_periods_manual*s%x(2)
        if(anomalous_diffusion_coefficient>0) call initialize_diffusion_reference(s)
        ! Measure the sign and frequency directly, before applying any OU kick.
        s%probe_x=s%x; s%probe_vp=s%vp; s%probe_vperp=s%vperp
        s%probe_tetr=s%tetr; s%probe_iface=s%iface; s%probe_initialized=.true.
        s%probe_dt=1e-5_dp/bg0(6)
        call orbit_timestep_gorilla(s%probe_x,s%probe_vp,s%probe_vperp,s%probe_dt, &
            s%probe_initialized,s%probe_tetr,s%probe_iface)
        if(s%probe_tetr==-1) error stop 'phase probe lost'
        s%phase_probe=mmode*theta_sfl(s%probe_x) &
            +nmode*n_field_periods_manual*s%probe_x(2)-s%alpha0
        s%phase_probe=atan2(sin(s%phase_probe),cos(s%phase_probe))/s%probe_dt
        s%predict_phase=bg0(4)+bg0(3)*s%v0
        s%phase_error=max(s%phase_error,abs(s%phase_probe-s%predict_phase)/bg0(6))
        if(im==1) then
            s%phase_row=[r0,s%v0/bg0(7),s%theta0,s%phase_probe/bg0(6),s%predict_phase/bg0(6)]
            s%probe_x=s%x; s%probe_vp=0.0_dp; s%probe_vperp=s%vperp
            s%probe_tetr=s%tetr; s%probe_iface=s%iface; s%probe_initialized=.true.
            call orbit_timestep_gorilla(s%probe_x,s%probe_vp,s%probe_vperp,s%probe_dt, &
                s%probe_initialized,s%probe_tetr,s%probe_iface)
            if(s%probe_tetr==-1) error stop 'zero-velocity phase probe lost'
            s%phase_probe=mmode*theta_sfl(s%probe_x)+nmode*n_field_periods_manual*s%probe_x(2)-s%alpha0
            s%phase_probe=atan2(sin(s%phase_probe),cos(s%phase_probe))/s%probe_dt
            s%zero_v_phase_error=max(s%zero_v_phase_error,abs(s%phase_probe-bg0(4))/bg0(6))
        end if
        if(ir==1) then
            s%probe_x=s%x; s%probe_vp=s%vp; s%probe_vperp=s%vperp
            s%probe_tetr=s%tetr; s%probe_iface=s%iface; s%probe_initialized=.true.
            s%probe_mu=-0.5_dp*s%vperp**2/bmod_func(s%x-tetra_physics(s%tetr)%x1,s%tetr)
            s%energy0=energy_tot_func([s%x-tetra_physics(s%tetr)%x1,s%vp],s%probe_mu,s%tetr)
            call orbit_timestep_gorilla(s%probe_x,s%probe_vp,s%probe_vperp,0.1_dp/bg0(6), &
                s%probe_initialized,s%probe_tetr,s%probe_iface)
            if(s%probe_tetr==-1) error stop 'conservation probe lost'
            s%energy1=energy_tot_func([s%probe_x-tetra_physics(s%probe_tetr)%x1,s%probe_vp],s%probe_mu,s%probe_tetr)
            s%max_energy_error=max(s%max_energy_error,abs(s%energy1-s%energy0)/abs(s%energy0))
            call orbit_timestep_gorilla(s%probe_x,s%probe_vp,s%probe_vperp,-0.1_dp/bg0(6), &
                s%probe_initialized,s%probe_tetr,s%probe_iface)
            if(s%probe_tetr==-1) error stop 'reversal probe lost'
            s%delta_x=s%probe_x-s%x
            s%delta_x(2)=R0_analytic_circ*atan2(sin(s%delta_x(2)),cos(s%delta_x(2)))
            s%max_reversal_error=max(s%max_reversal_error,sqrt(sum(s%delta_x**2))/r0, &
                abs(s%probe_vp-s%vp)/bg0(7),abs(s%probe_vperp-s%vperp)/bg0(7))
        end if
        s%p0=(s%v0/bg0(7))**2-1.0_dp
        s%local=0
        s%previous=0
        s%vc_phase=0
        s%t=0
        call sources(s)
        s%previous=s%present_source
        do while(s%t<end_t)
            s%phase_rate=max(abs(bg0(4))+abs(bg0(3))*max(abs(s%vc),bg0(7)), &
                abs(s%bgt(4))+abs(s%bgt(3))*max(abs(s%vp),s%bgt(7)))
            s%dt=min(nu_dt/max(bg0(6),s%bgt(6)),phase_step/max(s%phase_rate,1.0_dp),end_t-s%t)
            s%vc_old=s%vc
            call advance_response_step(s,-0.5_dp*s%dt,s%t)
            call random_number(s%rand2)
            s%zeta=sqrt(-2*log(max(s%rand2(1),tiny(1.0_dp))))*cos(2*pi*s%rand2(2))
            call exact_ou_velocity(s%vc,bg0(7)**2,bg0(6)*s%dt,s%zeta)
            call exact_ou_velocity(s%vp,s%bgt(7)**2,s%bgt(6)*s%dt,s%zeta)
            if(anomalous_diffusion_coefficient>0.0_dp) call diffuse_marker(s,s%dt)
            call advance_response_step(s,-0.5_dp*s%dt,s%t+0.5_dp*s%dt)
            s%vc_phase=s%vc_phase-s%dt*(bg0(4)+bg0(3)*(s%vc_old+s%vc)/2)
            s%t=s%t+s%dt
            call sources(s)
            s%local(:,1:2)=s%local(:,1:2)+0.5_dp*s%dt*(s%previous(:,1:2)+s%present_source(:,1:2))
            if(model==1) s%local(:,7)=s%local(:,7)+0.5_dp*s%dt*(s%previous(:,3)+s%present_source(:,3))
            s%previous=s%present_source
        end do
        if(s%n_respawn_used>0) s%n_respawn_markers=s%n_respawn_markers+1
        s%max_respawns_per_marker=max(s%max_respawns_per_marker,s%n_respawn_used)
        s%local(:,3:4)=s%local(:,1:2)*s%v0
        if(model==1) then
            s%local(:,5:6)=s%local(:,1:2)*s%p0
            s%local(:,8)=s%local(:,7)*s%v0
            s%local(:,9)=s%local(:,7)*s%p0
        end if
        result=s%local
        report=s%marker_report
    end subroutine
    subroutine reduce_report(report)
        type(marker_report),intent(in) :: report
        max_radial_excursion=max(max_radial_excursion,report%max_radial_excursion)
        phase_error=max(phase_error,report%phase_error)
        spawn_error=max(spawn_error,report%spawn_error)
        max_energy_error=max(max_energy_error,report%max_energy_error)
        max_reversal_error=max(max_reversal_error,report%max_reversal_error)
        zero_v_phase_error=max(zero_v_phase_error,report%zero_v_phase_error)
        max_respawns_per_marker=max(max_respawns_per_marker,report%max_respawns_per_marker)
        n_respawn_total=n_respawn_total+report%n_respawn_total
        n_respawn_markers=n_respawn_markers+report%n_respawn_markers
        n_respawn_inside=n_respawn_inside+report%n_respawn_inside
        n_respawn_outside=n_respawn_outside+report%n_respawn_outside
        diffusion_kicks=diffusion_kicks+report%diffusion_kicks
        diffusion_reflections=diffusion_reflections+report%diffusion_reflections
        diffusion_respawns=diffusion_respawns+report%diffusion_respawns
    end subroutine
    subroutine diffuse_marker(s,lag_step)
        ! Backward characteristics reverse Hamiltonian advection, not diffusion:
        ! the elliptic generator acts over positive lag. Reuse the applet's
        ! perpendicular tensor, geometric Ito drift and mesh traversal.
        type(marker_state),intent(inout) :: s
        real(dp),intent(in) :: lag_step
        integer :: count_steps,kick
        real(dp) :: delta_t,noise(3),reference_increment(3)
        logical :: reflected
        count_steps=max(1,ceiling(2*anomalous_diffusion_coefficient*lag_step/diffusion_step_cm**2))
        delta_t=lag_step/count_steps
        do kick=1,count_steps
            if(transport_boundary==1) then
                call anomalous_transport_displacement(s%x,s%tetr,s%iface,delta_t,s%vp,s%vperp, &
                    anomalous_diffusion_coefficient,rho_bounds=diffusion_bounds,recover_lost=.false., &
                    reflected=reflected,random_vector=noise)
                if(s%tetr==-1) error stop 'reflected diffusion displacement lost; no hidden axis recovery'
            else
                call anomalous_transport_displacement(s%x,s%tetr,s%iface,delta_t,s%vp,s%vperp, &
                    anomalous_diffusion_coefficient,recover_lost=.false.,random_vector=noise)
                reflected=.false.
                if(s%tetr==-1) then
                    ! This split kick consumes delta_t. On mesh exit discard
                    ! its untraced spatial remainder and apply the SAME reset
                    ! as a lost Hamiltonian path; subsequent kicks continue.
                    call respawn_response_marker(s,s%t+kick*delta_t,0.0_dp)
                    s%diffusion_respawns=s%diffusion_respawns+1
                end if
            end if
            if(.not.all(ieee_is_finite(s%x))) error stop 'invalid diffusion displacement'
            s%diffusion_kicks=s%diffusion_kicks+1
            if(reflected) s%diffusion_reflections=s%diffusion_reflections+1
            reference_increment=sqrt(2*delta_t)*matmul(s%reference_factor,noise)
            s%reference_radial_displacement=s%reference_radial_displacement &
                +dot_product(s%reference_radial_covector,reference_increment)
            s%reference_diffusion_phase=s%reference_diffusion_phase &
                +dot_product(s%reference_phase_covector,reference_increment)
        end do
        ! Mesh and magnetic moment are refreshed by the next polynomial push.
        ! Keep output labels, velocities, source integrals and coherent weights.
    end subroutine
    subroutine initialize_diffusion_reference(s)
        ! Flat cylindrical tangent frame. Recover its unit B direction from
        ! the exported parallel/perpendicular helical wave projections:
        ! kp=h_theta*(m/r)+h_phi*(n/Rref), ks=h_phi*(m/r)-h_theta*(n/Rref).
        ! This virtual reference has no radial boundary or geometry drift;
        ! the real reflected path remains in the traced difference.
        type(marker_state),intent(inout) :: s
        real(dp) :: a,c,htheta,hphi,norm,Rlocal,h(3),sintheta,costheta
        a=mmode/r0;c=nmode*n_field_periods_manual/R0_analytic_circ
        htheta=(a*bg0(3)-c*bg0(2))/(a*a+c*c)
        hphi=(c*bg0(3)+a*bg0(2))/(a*a+c*c)
        norm=sqrt(htheta*htheta+hphi*hphi);htheta=htheta/norm;hphi=hphi/norm
        sintheta=sin(s%theta0);costheta=cos(s%theta0);Rlocal=s%x(1)
        h=[-htheta*sintheta,hphi/Rlocal,htheta*costheta]
        call compute_diffusion_cholesky(h,Rlocal,anomalous_diffusion_coefficient,s%reference_factor)
        s%reference_radial_covector=[costheta,0.0_dp,sintheta]
        s%reference_phase_covector=bg0(2)*[-hphi*sintheta,-htheta*Rlocal,hphi*costheta]
    end subroutine
    subroutine advance_response_step(s,step,lag_start)
        ! Continue only the untraced Hamiltonian time after the common reset.
        type(marker_state),intent(inout) :: s
        real(dp),intent(in) :: step,lag_start
        real(dp) :: pending,remaining
        if(cylinder_only) return
        if(n_respawn_max==0) then
            call orbit_timestep_gorilla(s%x,s%vp,s%vperp,step,s%initialized,s%tetr,s%iface)
            if(s%tetr==-1) error stop 'orbit lost; respawn disabled'
            return
        end if
        pending=step
        do
            call orbit_timestep_gorilla(s%x,s%vp,s%vperp,pending,s%initialized,s%tetr,s%iface,remaining)
            if(s%tetr/=-1) return
            if(.not.ieee_is_finite(remaining).or..not.ieee_is_finite(s%vp).or. &
                .not.ieee_is_finite(s%vperp)) error stop 'invalid lost-particle state'
            if(remaining*step<0.0_dp.or.abs(remaining)>abs(pending)) &
                error stop 'invalid untraced pusher time'
            call respawn_response_marker(s,lag_start+abs(step-remaining),remaining)
            pending=remaining
            if(pending==0.0_dp) return
        end do
    end subroutine
    subroutine respawn_response_marker(s,lag,pending_time)
        ! Existing benchmark reset, shared by Hamiltonian and diffusion exits:
        ! redraw position uniformly in poloidal area in the unpadded KIM window
        ! and toroidal wedge. Keep velocities, elapsed lag, output labels,
        ! alpha0, virtual reference and all accumulated coherent channels.
        type(marker_state),intent(inout) :: s
        real(dp),intent(in) :: lag,pending_time
        real(dp),allocatable :: grown(:,:)
        real(dp) :: draw(3),old_radius,new_radius,old_vp,old_vperp
        real(dp) :: edge_flux,s_lo,s_hi,r_lo,r_hi
        integer :: attempt
        if(s%n_respawn_used>=n_respawn_max) error stop 'response respawn budget exhausted'
        if(.not.all(ieee_is_finite(s%x)).or..not.ieee_is_finite(s%vp).or. &
            .not.ieee_is_finite(s%vperp)) error stop 'invalid respawn state'
        r_lo=background(1,1); r_hi=r_lo+L
        ! Stable R0-sqrt(R0**2-r**2), in the benchmark's normalized s.
        edge_flux=a_analytic_circ**2/(R0_analytic_circ+sqrt(R0_analytic_circ**2-a_analytic_circ**2))
        s_lo=r_lo**2/(R0_analytic_circ+sqrt(R0_analytic_circ**2-r_lo**2))/edge_flux
        s_hi=r_hi**2/(R0_analytic_circ+sqrt(R0_analytic_circ**2-r_hi**2))/edge_flux
        old_radius=sqrt((s%x(1)-R0_analytic_circ)**2+s%x(3)**2)
        old_vp=s%vp; old_vperp=s%vperp
        do attempt=1,1000
            call random_number(draw)
            call draw_annulus_rphiz_analytic(s_lo,s_hi,draw,s%x)
            ! find_tetra, like the benchmark, rejects a rare placement at
            ! an edge. Select the entry face for the backward direction.
            call find_tetra(s%x,s%vp,s%vperp,s%tetr,s%iface,-1)
            if(s%tetr/=-1) exit
        end do
        if(s%tetr==-1) error stop 'could not respawn in response window'
        s%initialized=.true.
        s%n_respawn_used=s%n_respawn_used+1
        s%n_respawn_total=s%n_respawn_total+1
        if(old_radius>=a_analytic_circ*sqrt(sfc_s_min).and. &
            old_radius<=a_analytic_circ*sqrt(sfc_s_max)) then
            s%n_respawn_inside=s%n_respawn_inside+1
        else
            s%n_respawn_outside=s%n_respawn_outside+1
        end if
        new_radius=sqrt((s%x(1)-R0_analytic_circ)**2+s%x(3)**2)
        if(.not.allocated(s%events)) allocate(s%events(8,min(8,n_respawn_max)))
        if(s%n_respawn_total>size(s%events,2)) then
            allocate(grown(8,min(2*size(s%events,2),n_respawn_max)))
            grown(:,1:size(s%events,2))=s%events
            call move_alloc(grown,s%events)
        end if
        s%events(:,s%n_respawn_total)=[lag,pending_time, &
            old_radius,new_radius,old_vp,old_vperp,s%vp,s%vperp]
    end subroutine
    real(dp) function flux_radius(position,cell) result(radius)
        ! Match backgrounds to the unperturbed piecewise-linear poloidal flux,
        ! as in the existing applet, rather than to geometric distance from axis.
        real(dp),intent(in) :: position(3)
        integer,intent(in) :: cell
        real(dp) :: psi,f
        integer :: lo,hi,mid
        psi=abs(tetra_physics(cell)%Aphi1+sum(tetra_physics(cell)%gAphi*(position-tetra_physics(cell)%x1)))
        if(psi<flux_psi(1).or.psi>flux_psi(flux_points)) error stop 'flux radius outside mapping'
        lo=1; hi=flux_points
        do while(hi-lo>1)
            mid=(lo+hi)/2
            if(flux_psi(mid)>psi) then
                hi=mid
            else
                lo=mid
            end if
        end do
        f=(psi-flux_psi(lo))/(flux_psi(hi)-flux_psi(lo))
        radius=flux_r(lo)+f*(flux_r(hi)-flux_r(lo))
    end function
    real(dp) function theta_sfl(position) result(theta)
        ! Same exact analytic-circular SFL angle used by the core field module;
        ! evaluating it at x avoids an independent phase interpolation error.
        real(dp),intent(in) :: position(3)
        real(dp) :: geometric,eps,radius
        radius=sqrt((position(1)-R0_analytic_circ)**2+position(3)**2)
        geometric=atan2(position(3),position(1)-R0_analytic_circ)
        eps=radius/R0_analytic_circ
        theta=2*atan2(sqrt(1-eps)*sin(geometric/2),sqrt(1+eps)*cos(geometric/2))
    end function
    subroutine sample_background(r,bg)
        real(dp),intent(in) :: r
        real(dp),intent(out) :: bg(13)
        real(dp) :: index,f
        integer :: lo,hi
        index=modulo(r-background(1,1),L)*N/L
        lo=floor(index)+1
        hi=mod(lo,N)+1
        f=index-floor(index)
        bg=(1-f)*background(lo,:)+f*background(hi,:)
        bg(1)=r
    end subroutine
    subroutine sources(s)
        type(marker_state),intent(inout) :: s
        s%present_source(:,3)=0.0_dp
        if(cylinder_only) then
            s%rnow=r0
            s%bgt=bg0
            s%hphase=exp(ii*s%vc_phase)
        else
            s%rnow=flux_radius(s%x,s%tetr)
            call sample_background(s%rnow,s%bgt)
            s%alpha=mmode*theta_sfl(s%x)+nmode*n_field_periods_manual*s%x(2)
            s%hphase=exp(ii*(s%alpha-s%alpha0))
        end if
        s%max_radial_excursion=max(s%max_radial_excursion,abs(s%rnow-r0))
        s%coef=s%bgt(9)+s%bgt(10)*(1+0.5_dp*(s%vp/s%bgt(7))**2)
        s%source_phase=exp(ii*wave*s%rnow)*s%hphase
        s%present_source(:,1)=ii*clight*s%bgt(2)/s%bgt(13)*s%coef*s%source_phase
        s%present_source(:,2)=-s%vp/s%bgt(13)*s%coef*s%source_phase
        if(model==1) s%present_source(:,3)=s%bgt(6)*((s%vp/s%bgt(7))**2-1.0_dp)*s%source_phase
        s%cphase=exp(ii*s%vc_phase)
        if(anomalous_diffusion_coefficient>0) s%cphase=exp(ii*(s%vc_phase+s%reference_diffusion_phase))
        s%coef=coef0+0.5_dp*bg0(10)*(s%vc/bg0(7))**2
        s%source_phase=exp(ii*wave*r0)*s%cphase
        if(anomalous_diffusion_coefficient>0) &
            s%source_phase=exp(ii*wave*(r0+s%reference_radial_displacement))*s%cphase
        s%present_source(:,1)=s%present_source(:,1)-ii*clight*bg0(2)/bg0(13)*s%coef*s%source_phase
        s%present_source(:,2)=s%present_source(:,2)+s%vc/bg0(13)*s%coef*s%source_phase
        if(model==1) s%present_source(:,3)=s%present_source(:,3) &
            -bg0(6)*((s%vc/bg0(7))**2-1.0_dp)*s%source_phase
    end subroutine
    subroutine write_response(path,result)
        character(*),intent(in) :: path
        complex(dp),intent(in) :: result(:,:,:)
        integer :: u,i,jk,ic,irw
        open(newunit=u,file=path,status='new')
        if(size(result,3)==9) then
            write(u,'(a)') 'GK_ENERGY_BARE_V1'
        else
            write(u,'(a)') 'GK_RESPONSE_V1'
        end if
        write(u,*) M,N,model,mmode,nmode
        write(u,'(3es26.17e3)') L,rm,clight
        do i=1,N
            write(u,'(13es26.17e3)') background(i,:)
        end do
        do jk=1,size(result,3)
            do ic=1,d
                do irw=1,d
                    write(u,'(2es26.17e3)') real(result(irw,ic,jk)),aimag(result(irw,ic,jk))
                end do
            end do
        end do
        close(u)
    end subroutine
end program
