module orbit_ion_response_driver_mod
    ! Large-aspect-ratio ion response on actual equilibrium tetrahedral paths.
    ! Energy is closed at each GLOBAL magnetic moment before mu integration.
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use omp_lib, only: omp_get_wtime
    use cylindrical_ion_mu_response_mod, only: ion_mu_source,ion_mu_deposit,assemble_cylindrical_ion_mu, &
        assemble_diffused_ion_mu
    use cylindrical_ion_response_mod, only: perpendicular_rule
    use cylindrical_electron_response_mod, only: close_energy_response
    use orbit_timestep_gorilla_mod, only: initialize_gorilla,orbit_timestep_gorilla
    use gorilla_settings_mod, only: load_gorilla_inp,ipusher,coord_system, &
        boole_helical_pert,boole_pert_from_mephit,boole_axi_noise_vector_pot, &
        boole_axi_noise_elec_pot,boole_non_axi_noise_vector_pot,boole_strong_electric_field,boole_magnetic_drift
    use tetra_grid_settings_mod, only: load_tetra_grid_inp,R0_analytic_circ,n_field_periods_manual,grid_kind, &
        a_analytic_circ,sfc_s_min,sfc_s_max,n3
    use tetra_physics_mod, only: tetra_physics,particle_mass,particle_charge
    use supporting_functions_mod, only: bmod_func,energy_tot_func
    use field_analytic_circ_mod, only: psi_pol_analytic_circ
    use collis_ions, only: exact_ou_velocity,getran
    use pusher_tetra_func_mod, only: pusher_handover2neighbour
    use anomalous_transport_displacement_mod, only: anomalous_transport_displacement,compute_diffusion_cholesky
    implicit none
    private
    public :: run_orbit_ion_response
    real(dp),parameter :: pi=acos(-1.0_dp)
    complex(dp),parameter :: ii=(0.0_dp,1.0_dp)
    character(1024) :: input_path,output_path,tag,arg
    integer :: M,N,model,mmode,nmode,species,unit,j,k,d,imu,ir,im,ib,col,refinement=1
    integer :: radial_samples=16,markers=32,batches=4,n_mu=8,seed=18371
    real(dp) :: L,rm,clight,lag_nu=8,nu_dt=.05_dp,phase_step=.15_dp,mu_reference_scale=1,prefix_lag_nu=0
    ! Uniform shift of all injection radii [cm]; zero keeps the legacy sample positions.
    real(dp) :: radial_sample_shift_cm=0
    logical :: cylinder_only=.false.
    ! Constant perpendicular anomalous diffusion [cm^2/s] of the guiding centre,
    ! split into kicks of nominal RMS length diffusion_step_cm. Zero is legacy.
    real(dp) :: anomalous_diffusion_coefficient=0,diffusion_step_cm=.1_dp
    ! Cylinder-only check: deposit the diffused virtual path itself instead of
    ! the paired difference, to compare its raw channels with the reference.
    logical :: diffusion_reference_check=.false.
    ! Paired estimator: trace every marker with and without diffusion on shared
    ! Hamiltonian/collision draws (separate deterministic stream for the kicks),
    ! writing the D=0 partner to *_D0* files for paired differences.
    logical :: paired_diffusion=.false.
    ! Contiguous radial-sample subset for parallel processes; 0 means all.
    ! Raw per-mu corrections (ion_paired_mu.dat) of disjoint subsets add exactly.
    integer :: radial_sample_first=1,radial_sample_last=0
    ! Reseed the intrinsic stream per (mu node, radial sample) so that disjoint
    ! subsets reproduce a full run; required for subsets, default legacy stream.
    logical :: sample_streams=.false.
    ! Write the raw per-mu channels as unformatted stream (ion_paired_mu.bin)
    ! instead of text; needed for large Fourier bases.
    logical :: raw_binary=.false.
    integer :: raw_unit
    integer :: diffusion_relocations=0
    ! Diagnostic trace of one marker (mu node, radial sample, marker): every
    ! kick and pusher step with signed face distances; 0 disables. Output only.
    integer :: debug_marker(3)=0,debug_unit=0
    namelist /ion_orbit_response/ radial_samples,markers,batches,n_mu,seed,lag_nu,nu_dt,phase_step, &
        mu_reference_scale,cylinder_only,prefix_lag_nu,radial_sample_shift_cm, &
        anomalous_diffusion_coefficient,diffusion_step_cm,diffusion_reference_check,paired_diffusion, &
        radial_sample_first,radial_sample_last,sample_streams,raw_binary,debug_marker
    complex(dp),allocatable :: exact0(:,:,:),bare0(:,:,:),cylinder_closed0(:,:,:)
    complex(dp),allocatable :: corrections0(:,:,:,:),mean0(:,:,:),batches0(:,:,:,:),integrated_exact0(:,:,:)
    complex(dp),allocatable :: prefix_corrections0(:,:,:,:),prefix_mean0(:,:,:),prefix_batches0(:,:,:,:)
    complex(dp),allocatable :: sample0(:,:,:),prefix_sample0(:,:,:)
    integer :: closure_unit0,prefix_closure_unit0,paired_unit0
    real(dp) :: diffusion_bounds(2)=0
    integer :: diffusion_kicks=0,diffusion_reflections=0
    real(dp),allocatable :: bg(:,:),wave(:),nodes(:),weights(:)
    integer,allocatable :: seed_array(:)
    real(dp) :: bg0(14),reference_bg(14),mu_ref,lambda,r0,closure_error,rcond,mesh_seconds,start_time
    real(dp) :: max_energy_error=0,max_reversal_error=0,max_mu_error=0,max_excursion=0,max_phase_difference=0
    real(dp) :: mass_from_background,charge_from_background
    complex(dp),allocatable :: exact(:,:,:),bare(:,:,:),cylinder_closed(:,:,:),closed(:,:,:)
    complex(dp),allocatable :: corrections(:,:,:,:),mean_result(:,:,:),batch_results(:,:,:,:)
    complex(dp),allocatable :: integrated_exact(:,:,:)
    complex(dp),allocatable :: projection(:),sample_result(:,:,:)
    integer,parameter :: flux_points=8193
    real(dp) :: flux_r(flux_points),flux_psi(flux_points)
    complex(dp),allocatable :: prefix_corrections(:,:,:,:),prefix_mean(:,:,:),prefix_batches(:,:,:,:)
    complex(dp),allocatable :: prefix_sample(:,:,:)
    integer :: progress,diagnostics,closure_unit,prefix_closure_unit,paired_unit,prefix_checkpoint_count=0
contains
    subroutine run_orbit_ion_response()
    call get_command_argument(1,input_path)
    call get_command_argument(2,output_path)
    call get_command_argument(3,arg)
    if(len_trim(input_path)==0.or.len_trim(output_path)==0) &
        error stop 'usage: gorilla_orbit_ion_response.x background response [refinement]'
    if(len_trim(arg)>0) read(arg,*) refinement
    open(newunit=unit,file=trim(input_path),status='old')
    read(unit,'(a)') tag
    if(trim(tag)/='GK_ION_BACKGROUND_V1') error stop 'wrong ion background format'
    read(unit,*) M,N,model,mmode,nmode,species
    read(unit,*) L,rm,clight
    if(M<0.or.N<2*M+1.or.L<=0.or.clight<=0.or.species<1.or.refinement<1) &
        error stop 'invalid ion orbit metadata'
    if(model/=0.and.model/=1) error stop 'ion orbit provider supports OU models 0/1'
    allocate(bg(N,14));do j=1,N;read(unit,*) bg(j,:);end do;close(unit)
    if(.not.all(ieee_is_finite(bg))) error stop 'nonfinite ion background'
    if(any(bg(:,5:7)<=0).or.any(bg(:,8)==0).or.any(bg(:,13)<=0).or.any(bg(:,14)<0)) &
        error stop 'invalid ion thermal parameters'
    open(newunit=unit,file='ion_orbit_response.nml',status='old')
    read(unit,nml=ion_orbit_response);close(unit)
    if(radial_samples<2*M+1.or.markers<1.or.batches<2.or.mod(markers,batches)/=0.or.n_mu<2) &
        error stop 'invalid ion orbit quadrature/batches'
    if(.not.all(ieee_is_finite([lag_nu,nu_dt,phase_step,mu_reference_scale,prefix_lag_nu,radial_sample_shift_cm]))) &
        error stop 'nonfinite ion orbit settings'
    if(2*abs(radial_sample_shift_cm)*radial_samples>=L) &
        error stop 'radial_sample_shift_cm must stay inside half a sample spacing'
    if(.not.ieee_is_finite(anomalous_diffusion_coefficient).or.anomalous_diffusion_coefficient<0) &
        error stop 'invalid anomalous diffusion coefficient'
    if(anomalous_diffusion_coefficient>0.and.(.not.ieee_is_finite(diffusion_step_cm).or.diffusion_step_cm<=0)) &
        error stop 'invalid diffusion substep length'
    if(diffusion_reference_check.and.(.not.cylinder_only.or.prefix_lag_nu<=0)) &
        error stop 'diffusion reference check requires cylinder_only and a prefix window'
    if(paired_diffusion.and.(anomalous_diffusion_coefficient<=0.or.diffusion_reference_check.or.prefix_lag_nu<=0)) &
        error stop 'paired diffusion requires D>0, a prefix window and no reference check'
    if(radial_sample_last==0) radial_sample_last=radial_samples
    if(radial_sample_first<1.or.radial_sample_last>radial_samples.or.radial_sample_first>radial_sample_last) &
        error stop 'invalid radial sample subset'
    if((radial_sample_first/=1.or.radial_sample_last/=radial_samples).and.(prefix_lag_nu<=0.or..not.sample_streams)) &
        error stop 'radial subsets require a prefix window (raw per-mu channels) and sample_streams'
    if(min(lag_nu,nu_dt,phase_step,mu_reference_scale)<=0) error stop 'invalid ion orbit step/lag'
    if(prefix_lag_nu<0.or.prefix_lag_nu>lag_nu) error stop 'prefix must be zero or within total lag'
    call sample_background(rm,reference_bg)
    mu_ref=mu_reference_scale*reference_bg(7)**2/reference_bg(13)
    call perpendicular_rule(n_mu,nodes,weights)
    call random_seed(size=j);allocate(seed_array(j));seed_array=seed
    call random_seed(put=seed_array)
    mesh_seconds=0
    if(.not.cylinder_only) then
        call load_gorilla_inp();call load_tetra_grid_inp()
        if(grid_kind/=5.or.coord_system/=1.or.ipusher/=2) &
            error stop 'ion benchmark requires analytic-circular cylindrical-coordinate polynomial pusher'
        if(boole_helical_pert.or.boole_pert_from_mephit.or.boole_axi_noise_vector_pot.or. &
            boole_axi_noise_elec_pot.or.boole_non_axi_noise_vector_pot) &
            error stop 'linear ion response requires an unperturbed equilibrium orbit field'
        if(.not.boole_strong_electric_field) error stop 'ion benchmark requires the matched equilibrium Er'
        start_time=omp_get_wtime()
        call initialize_gorilla(ipert_in=0)
        mesh_seconds=omp_get_wtime()-start_time
        do j=1,N
            mass_from_background=bg(j,12)*1.6022e-12_dp/bg(j,7)**2
            charge_from_background=bg(j,8)*mass_from_background*clight/bg(j,13)
            ! Arithmetic consistency of the actual exported mass/charge, not a physics tolerance.
            if(abs(mass_from_background/particle_mass-1)>256*epsilon(1.0_dp).or. &
                abs(charge_from_background/particle_charge-1)>256*epsilon(1.0_dp)) &
                error stop 'GORILLA pusher species differs from exported ion species'
        end do
        do j=1,flux_points
            flux_r(j)=1+49*real(j-1,dp)/(flux_points-1)
            flux_psi(j)=abs(psi_pol_analytic_circ(flux_r(j)))
        end do
        if(anomalous_diffusion_coefficient>0) then
            ! Reflect inside the mesh-contained annulus (same construction as the
            ! electron provider): inner vertex circle and outer polygon apothem.
            diffusion_bounds=sqrt([sfc_s_min,sfc_s_max]*(R0_analytic_circ-sqrt(R0_analytic_circ**2 &
                -a_analytic_circ**2))*(2*R0_analytic_circ-[sfc_s_min,sfc_s_max]*(R0_analytic_circ &
                -sqrt(R0_analytic_circ**2-a_analytic_circ**2))))
            diffusion_bounds(2)=diffusion_bounds(2)*cos(pi/n3)
            if(diffusion_bounds(1)>=diffusion_bounds(2)) error stop 'mesh-contained reflecting annulus is empty'
        end if
    end if
    d=2*M+1
    allocate(wave(d),corrections(d,d,9,batches),mean_result(d,d,4),batch_results(d,d,4,batches))
    allocate(projection(d),sample_result(d,d,9),integrated_exact(d,d,4))
    wave=2*pi*real([(j-M-1,j=1,d)],dp)/L
    mean_result=0;batch_results=0;integrated_exact=0
    allocate(prefix_sample(d,d,9))
    if(prefix_lag_nu>0) then
        allocate(prefix_corrections(d,d,9,batches),prefix_mean(d,d,4),prefix_batches(d,d,4,batches))
        prefix_mean=0;prefix_batches=0
        open(newunit=prefix_closure_unit,file='ion_prefix_mu_closure.dat',status='new')
        write(prefix_closure_unit,'(a)') '# node lambda weight mean_closure_residual mean_rcond'
        if(raw_binary) then
            ! Layout per mu node: int32 node,d,batches; complex128 bare(d,d,9),
            ! exact(d,d,4), prefix(d,d,9,batches), total(d,d,9,batches).
            open(newunit=raw_unit,file='ion_paired_mu.bin',status='new',access='stream',form='unformatted')
        else
            open(newunit=paired_unit,file='ion_paired_mu.dat',status='new')
            write(paired_unit,'(a)') '# node batch window channel column row real imag'
            write(paired_unit,'(a)') '# window 0=bare, 1=prefix correction, 2=total correction, 3=exact closed'
        end if
    end if
    if(paired_diffusion) then
        allocate(corrections0(d,d,9,batches),mean0(d,d,4),batches0(d,d,4,batches),integrated_exact0(d,d,4))
        allocate(prefix_corrections0(d,d,9,batches),prefix_mean0(d,d,4),prefix_batches0(d,d,4,batches))
        allocate(sample0(d,d,9),prefix_sample0(d,d,9))
        mean0=0;batches0=0;integrated_exact0=0;prefix_mean0=0;prefix_batches0=0
        open(newunit=closure_unit0,file='ion_mu_closure_D0.dat',status='new')
        write(closure_unit0,'(a)') '# node lambda weight mean_closure_residual mean_rcond'
        open(newunit=prefix_closure_unit0,file='ion_prefix_mu_closure_D0.dat',status='new')
        write(prefix_closure_unit0,'(a)') '# node lambda weight mean_closure_residual mean_rcond'
        open(newunit=paired_unit0,file='ion_paired_mu_D0.dat',status='new')
        write(paired_unit0,'(a)') '# node batch window channel column row real imag (D=0 partner)'
    end if
    open(newunit=progress,file='ion_orbit_progress.dat',status='new')
    write(progress,'(a)') '# mu_node radial_sample elapsed_seconds'
    open(newunit=closure_unit,file='ion_mu_closure.dat',status='new')
    write(closure_unit,'(a)') '# node lambda weight mean_closure_residual mean_rcond'
    start_time=omp_get_wtime()
    do imu=1,n_mu
        lambda=nodes(imu)
        call assemble_cylindrical_ion_mu(bg,M,L,clight,lambda,mu_ref,model,exact,bare,refinement)
        if(anomalous_diffusion_coefficient>0) then
            ! Diffused reference with the same structure as D=0: locally restored
            ! exact part and bare channels for the paired field-space closure.
            if(paired_diffusion) then
                call move_alloc(exact,exact0);call move_alloc(bare,bare0)
                integrated_exact0=integrated_exact0+weights(imu)*exact0
                if(model==1) call close_energy_response(bare0,cylinder_closed0,closure_error,rcond)
                corrections0=0;prefix_corrections0=0
            end if
            if(allocated(exact)) deallocate(exact,bare)
            call assemble_diffused_ion_mu(bg,M,L,clight,lambda,mu_ref,model,anomalous_diffusion_coefficient,exact,bare)
        end if
        integrated_exact=integrated_exact+weights(imu)*exact
        if(model==1) call close_energy_response(bare,cylinder_closed,closure_error,rcond)
        corrections=0
        if(prefix_lag_nu>0) prefix_corrections=0
        do ir=radial_sample_first,radial_sample_last
            r0=bg(1,1)+(ir-.5_dp)*L/radial_samples+radial_sample_shift_cm
            if(sample_streams) then
                seed_array=seed+1000003*imu+7919*ir
                call random_seed(put=seed_array)
            end if
            call sample_background(r0,bg0)
            projection=exp(-ii*wave*r0)/real(radial_samples*markers/batches,dp)
            do im=1,markers
                if(paired_diffusion) then
                    call trace_marker_paired(ir,im,sample_result,prefix_sample,sample0,prefix_sample0)
                else
                    call trace_marker(ir,im,sample_result,prefix_sample)
                end if
                ib=mod(im-1,batches)+1
                if(paired_diffusion) then
                    do k=1,9
                        do col=1,d
                            corrections0(:,col,k,ib)=corrections0(:,col,k,ib)+projection*sample0(:,col,k)
                            prefix_corrections0(:,col,k,ib)=prefix_corrections0(:,col,k,ib) &
                                +projection*prefix_sample0(:,col,k)
                        end do
                    end do
                end if
                do k=1,9
                    do col=1,d
                        corrections(:,col,k,ib)=corrections(:,col,k,ib)+projection*sample_result(:,col,k)
                        if(prefix_lag_nu>0) prefix_corrections(:,col,k,ib)= &
                            prefix_corrections(:,col,k,ib)+projection*prefix_sample(:,col,k)
                    end do
                end do
            end do
            write(progress,'(2i8,es24.16)') imu,ir,omp_get_wtime()-start_time
            flush(progress)
        end do
        closure_error=0;rcond=1
        if(model==1) then
            call close_energy_response(bare+sum(corrections,dim=4)/batches,closed,closure_error,rcond)
            mean_result=mean_result+weights(imu)*(exact+closed-cylinder_closed)
            deallocate(closed)
        else
            mean_result=mean_result+weights(imu)*(exact+sum(corrections(:,:,1:4,:),dim=4)/batches)
        end if
        write(closure_unit,'(i8,4es24.16)') imu,lambda,weights(imu),closure_error,rcond
        flush(closure_unit)
        do ib=1,batches
            if(model==1) then
                call close_energy_response(bare+corrections(:,:,:,ib),closed,closure_error,rcond)
                batch_results(:,:,:,ib)=batch_results(:,:,:,ib)+weights(imu)*(exact+closed-cylinder_closed)
                deallocate(closed)
            else
                batch_results(:,:,:,ib)=batch_results(:,:,:,ib)+weights(imu)*(exact+corrections(:,:,1:4,ib))
            end if
        end do
        if(prefix_lag_nu>0) then
            closure_error=0;rcond=1
            if(model==1) then
                call close_energy_response(bare+sum(prefix_corrections,dim=4)/batches,closed,closure_error,rcond)
                prefix_mean=prefix_mean+weights(imu)*(exact+closed-cylinder_closed)
                deallocate(closed)
            else
                prefix_mean=prefix_mean+weights(imu)*(exact+sum(prefix_corrections(:,:,1:4,:),dim=4)/batches)
            end if
            write(prefix_closure_unit,'(i8,4es24.16)') imu,lambda,weights(imu),closure_error,rcond
            flush(prefix_closure_unit)
            if(raw_binary) then
                write(raw_unit) int(imu,4),int(d,4),int(batches,4),bare,exact,prefix_corrections,corrections
                flush(raw_unit)
            else
                call write_mu_channels(0,0,bare)
                call write_mu_channels(0,3,exact)
            end if
            do ib=1,batches
                if(.not.raw_binary) then
                    call write_mu_channels(ib,1,prefix_corrections(:,:,:,ib))
                    call write_mu_channels(ib,2,corrections(:,:,:,ib))
                end if
                if(model==1) then
                    call close_energy_response(bare+prefix_corrections(:,:,:,ib),closed,closure_error,rcond)
                    prefix_batches(:,:,:,ib)=prefix_batches(:,:,:,ib)+weights(imu)*(exact+closed-cylinder_closed)
                    deallocate(closed)
                else
                    prefix_batches(:,:,:,ib)=prefix_batches(:,:,:,ib)+weights(imu)*(exact+prefix_corrections(:,:,1:4,ib))
                end if
            end do
            if(.not.raw_binary) flush(paired_unit)
        end if
        if(paired_diffusion) call close_paired_partner()
        deallocate(exact,bare)
        if(allocated(cylinder_closed)) deallocate(cylinder_closed)
        print *, 'Completed ion magnetic-moment node ',imu,'/',n_mu
        flush(6)
    end do
    close(progress);close(closure_unit)
    call write_response(trim(output_path),mean_result)
    call write_response('ion_orbit_cylinder_reference.dat',integrated_exact)
    do ib=1,batches
        write(arg,'(a,i0,a)') 'ion_orbit_batch_',ib,'.dat'
        call write_response(trim(arg),batch_results(:,:,:,ib))
    end do
    if(paired_diffusion) then
        close(closure_unit0);close(prefix_closure_unit0);close(paired_unit0)
        call write_response('ion_orbit_D0_response.dat',mean0)
        call write_response('ion_orbit_D0_prefix_response.dat',prefix_mean0)
        call write_response('ion_orbit_D0_cylinder_reference.dat',integrated_exact0)
        do ib=1,batches
            write(arg,'(a,i0,a)') 'ion_orbit_D0_batch_',ib,'.dat'
            call write_response(trim(arg),batches0(:,:,:,ib))
            write(arg,'(a,i0,a)') 'ion_orbit_D0_prefix_batch_',ib,'.dat'
            call write_response(trim(arg),prefix_batches0(:,:,:,ib))
        end do
    end if
    if(prefix_lag_nu>0) then
        close(prefix_closure_unit)
        if(raw_binary) then
            close(raw_unit)
        else
            close(paired_unit)
        end if
        call write_response('ion_orbit_prefix_response.dat',prefix_mean)
        do ib=1,batches
            write(arg,'(a,i0,a)') 'ion_orbit_prefix_batch_',ib,'.dat'
            call write_response(trim(arg),prefix_batches(:,:,:,ib))
        end do
    end if
    open(newunit=diagnostics,file='ion_orbit_diagnostics.dat',status='new')
    write(diagnostics,'(a,es24.16)') 'mesh_initialization_seconds ',mesh_seconds
    write(diagnostics,'(a,es24.16)') 'response_seconds ',omp_get_wtime()-start_time
    write(diagnostics,'(a,es24.16)') 'mu_reference_cm2_s2_per_G ',mu_ref
    write(diagnostics,'(a,es24.16)') 'max_relative_magnetic_moment_error ',max_mu_error
    write(diagnostics,'(a,es24.16)') 'max_collisionless_energy_error ',max_energy_error
    write(diagnostics,'(a,es24.16)') 'max_collisionless_reversal_error ',max_reversal_error
    write(diagnostics,'(a,es24.16)') 'max_radial_excursion_cm ',max_excursion
    write(diagnostics,'(a,es24.16)') 'max_phase_frequency_difference_over_nu ',max_phase_difference
    ! Written only for the drift-free diagnostic, so default outputs are unchanged.
    if(.not.cylinder_only.and..not.boole_magnetic_drift) write(diagnostics,'(a,i12)') 'magnetic_drift_disabled ',1
    ! Written only with diffusion, so legacy outputs are unchanged.
    if(anomalous_diffusion_coefficient>0) then
        write(diagnostics,'(a,es24.16)') 'anomalous_diffusion_cm2_per_s ',anomalous_diffusion_coefficient
        write(diagnostics,'(a,es24.16)') 'diffusion_step_cm ',diffusion_step_cm
        write(diagnostics,'(a,i12)') 'diffusion_kicks ',diffusion_kicks
        write(diagnostics,'(a,i12)') 'diffusion_reflections ',diffusion_reflections
        if(.not.cylinder_only) then
            write(diagnostics,'(a,es24.16)') 'diffusion_lower_boundary_cm ',diffusion_bounds(1)
            write(diagnostics,'(a,es24.16)') 'diffusion_upper_boundary_cm ',diffusion_bounds(2)
        end if
    end if
    if(diffusion_reference_check) write(diagnostics,'(a,i12)') 'diffusion_reference_check ',1
    if(paired_diffusion) write(diagnostics,'(a,i12)') 'paired_diffusion ',1
    if(sample_streams) write(diagnostics,'(a,i12)') 'sample_streams ',1
    if(anomalous_diffusion_coefficient>0) write(diagnostics,'(a,i12)') 'diffusion_relocations ',diffusion_relocations
    if(radial_sample_first/=1.or.radial_sample_last/=radial_samples) then
        write(diagnostics,'(a,i12)') 'radial_sample_first ',radial_sample_first
        write(diagnostics,'(a,i12)') 'radial_sample_last ',radial_sample_last
    end if
    if(prefix_lag_nu>0) then
        if(prefix_checkpoint_count/=n_mu*(radial_sample_last-radial_sample_first+1)*markers) &
            error stop 'missing marker prefix checkpoint'
        write(diagnostics,'(a,i12)') 'prefix_checkpoint_count ',prefix_checkpoint_count
    end if
    close(diagnostics)
    contains
        subroutine close_paired_partner()
            ! Same per-mu closure and outputs as the main estimator, for the D=0 partner.
            complex(dp),allocatable :: partner(:,:,:)
            integer :: b
            closure_error=0;rcond=1
            if(model==1) then
                call close_energy_response(bare0+sum(corrections0,dim=4)/batches,partner,closure_error,rcond)
                mean0=mean0+weights(imu)*(exact0+partner-cylinder_closed0);deallocate(partner)
            else
                mean0=mean0+weights(imu)*(exact0+sum(corrections0(:,:,1:4,:),dim=4)/batches)
            end if
            write(closure_unit0,'(i8,4es24.16)') imu,lambda,weights(imu),closure_error,rcond
            closure_error=0;rcond=1
            if(model==1) then
                call close_energy_response(bare0+sum(prefix_corrections0,dim=4)/batches,partner,closure_error,rcond)
                prefix_mean0=prefix_mean0+weights(imu)*(exact0+partner-cylinder_closed0);deallocate(partner)
            else
                prefix_mean0=prefix_mean0+weights(imu)*(exact0+sum(prefix_corrections0(:,:,1:4,:),dim=4)/batches)
            end if
            write(prefix_closure_unit0,'(i8,4es24.16)') imu,lambda,weights(imu),closure_error,rcond
            call write_partner_channels(0,0,bare0);call write_partner_channels(0,3,exact0)
            do b=1,batches
                call write_partner_channels(b,1,prefix_corrections0(:,:,:,b))
                call write_partner_channels(b,2,corrections0(:,:,:,b))
                if(model==1) then
                    call close_energy_response(bare0+corrections0(:,:,:,b),partner,closure_error,rcond)
                    batches0(:,:,:,b)=batches0(:,:,:,b)+weights(imu)*(exact0+partner-cylinder_closed0)
                    deallocate(partner)
                    call close_energy_response(bare0+prefix_corrections0(:,:,:,b),partner,closure_error,rcond)
                    prefix_batches0(:,:,:,b)=prefix_batches0(:,:,:,b)+weights(imu)*(exact0+partner-cylinder_closed0)
                    deallocate(partner)
                else
                    batches0(:,:,:,b)=batches0(:,:,:,b)+weights(imu)*(exact0+corrections0(:,:,1:4,b))
                    prefix_batches0(:,:,:,b)=prefix_batches0(:,:,:,b)+weights(imu)*(exact0+prefix_corrections0(:,:,1:4,b))
                end if
            end do
            flush(closure_unit0);flush(prefix_closure_unit0);flush(paired_unit0)
            deallocate(exact0,bare0)
            if(allocated(cylinder_closed0)) deallocate(cylinder_closed0)
        end subroutine
        subroutine write_partner_channels(batch_index,window_index,values)
            integer,intent(in) :: batch_index,window_index
            complex(dp),intent(in) :: values(:,:,:)
            integer :: row,column,channel
            do channel=1,size(values,3)
                do column=1,d
                    do row=1,d
                        write(paired_unit0,'(6i6,2es26.17e3)') imu,batch_index,window_index,channel,column,row, &
                            real(values(row,column,channel)),aimag(values(row,column,channel))
                    end do
                end do
            end do
        end subroutine
    end subroutine run_orbit_ion_response
    subroutine write_mu_channels(batch_index,window_index,values)
        integer,intent(in) :: batch_index,window_index
        complex(dp),intent(in) :: values(:,:,:)
        integer :: row,column,channel
        do channel=1,size(values,3)
            do column=1,d
                do row=1,d
                    write(paired_unit,'(6i6,2es26.17e3)') imu,batch_index,window_index,channel,column,row, &
                        real(values(row,column,channel)),aimag(values(row,column,channel))
                end do
            end do
        end do
    end subroutine write_mu_channels
    subroutine trace_marker(radial_index,marker_index,result,prefix_result)
        integer,intent(in) :: radial_index,marker_index
        complex(dp),intent(out) :: result(d,d,9),prefix_result(d,d,9)
        real(dp) :: x(3),vpar,vperp,v0,vc,vc_old,t,dt,rand2(2),theta0,alpha0,alpha,phase_c
        real(dp) :: current_bg(14),rnow,geometric,field,field0,mu,phase_rate,end_time,zeta
        real(dp) :: prefix_time,target_time
        logical :: prefix_saved
        real(dp) :: probe_x(3),probe_vpar,probe_vperp,energy0,energy1,probe_dt,phase0,delta(3)
        complex(dp) :: actual_source(d,3),ref_source(d,3),old_actual(d,3),old_ref(d,3)
        complex(dp) :: real_integral(d,3),ref_integral(d,3),real_deposit(d,d,9),ref_deposit(d,d,9)
        logical :: initialized,probe_initialized
        integer :: cell,face,iteration,probe_cell,probe_face
        real(dp) :: ref_dr,ref_dphase,ref_radial(3),ref_helical(3)
        logical :: debugging
        debugging=all([imu,radial_index,marker_index]==debug_marker)
        if(debugging.and.debug_unit==0) open(newunit=debug_unit,file='ion_debug_marker.dat',status='new')
        call random_number(rand2)
        v0=bg0(7)*sqrt(-2*log(max(rand2(1),tiny(1.0_dp))))*cos(2*pi*rand2(2))
        vpar=v0;vc=v0;t=0;phase_c=0;mu=lambda*mu_ref
        ref_dr=0;ref_dphase=0
        rnow=r0;current_bg=bg0;field0=bg0(13)
        if(.not.cylinder_only) then
            call random_number(rand2);theta0=2*pi*rand2(1);geometric=r0
            vperp=sqrt(2*mu*bg0(13))
            do iteration=1,5
                x=[R0_analytic_circ+geometric*cos(theta0), &
                    2*pi*rand2(2)/n_field_periods_manual,geometric*sin(theta0)]
                initialized=.false.;cell=0;face=0
                call orbit_timestep_gorilla(x,vpar,vperp,0.0_dp,initialized,cell,face)
                if(cell==-1) error stop 'ion injection outside mesh'
                rnow=flux_radius(x,cell)
                if(abs(rnow-r0)<1e-9_dp) exit
                geometric=geometric+r0-rnow
            end do
            if(abs(rnow-r0)>1e-8_dp) error stop 'ion injection flux label did not converge'
            field0=bmod_func(x-tetra_physics(cell)%x1,cell)
            vperp=sqrt(2*mu*field0)
            initialized=.false.
            call orbit_timestep_gorilla(x,vpar,vperp,0.0_dp,initialized,cell,face)
            alpha0=mmode*theta_sfl(x)+nmode*n_field_periods_manual*x(2)
            if(anomalous_diffusion_coefficient>0) call set_reference_frame(theta0,x(1))
            if(radial_index==1.and.marker_index==1) then
                probe_x=x;probe_vpar=vpar;probe_vperp=vperp
                probe_cell=cell;probe_face=face;probe_initialized=.true.
                energy0=energy_tot_func([x-tetra_physics(cell)%x1,vpar],-mu,cell)
                probe_dt=.01_dp/bg0(6)
                call orbit_timestep_gorilla(probe_x,probe_vpar,probe_vperp,probe_dt, &
                    probe_initialized,probe_cell,probe_face)
                if(probe_cell==-1) error stop 'ion conservation probe lost'
                energy1=energy_tot_func([probe_x-tetra_physics(probe_cell)%x1,probe_vpar],-mu,probe_cell)
                max_energy_error=max(max_energy_error,abs(energy1-energy0)/max(abs(energy0),tiny(1.0_dp)))
                phase0=mmode*theta_sfl(probe_x)+nmode*n_field_periods_manual*probe_x(2)-alpha0
                max_phase_difference=max(max_phase_difference, &
                    abs(atan2(sin(phase0),cos(phase0))/probe_dt-bg0(4)-bg0(3)*v0)/bg0(6))
                call orbit_timestep_gorilla(probe_x,probe_vpar,probe_vperp,-probe_dt, &
                    probe_initialized,probe_cell,probe_face)
                if(probe_cell==-1) error stop 'ion reversal probe lost'
                delta=probe_x-x
                delta(2)=R0_analytic_circ/n_field_periods_manual &
                    *atan2(sin(n_field_periods_manual*delta(2)),cos(n_field_periods_manual*delta(2)))
                max_reversal_error=max(max_reversal_error,sqrt(sum(delta**2))/r0, &
                    abs(probe_vpar-vpar)/bg0(7),abs(probe_vperp-vperp)/bg0(7))
            end if
        end if
        real_integral=0;ref_integral=0
        call sources();old_actual=actual_source;old_ref=ref_source
        end_time=lag_nu/bg0(6)
        prefix_time=prefix_lag_nu/bg0(6)
        prefix_saved=prefix_lag_nu<=0
        prefix_result=0
        do while(t<end_time)
            phase_rate=max(abs(bg0(4))+abs(bg0(3))*max(abs(vc),bg0(7)), &
                abs(current_bg(4))+abs(current_bg(3))*max(abs(vpar),current_bg(7)))
            target_time=end_time
            if(.not.prefix_saved) target_time=prefix_time
            dt=min(nu_dt/max(bg0(6),current_bg(6)),phase_step/max(phase_rate,1.0_dp),target_time-t)
            vc_old=vc
            call advance(-dt/2)
            call random_number(rand2)
            zeta=sqrt(-2*log(max(rand2(1),tiny(1.0_dp))))*cos(2*pi*rand2(2))
            call exact_ou_velocity(vc,bg0(7)**2,bg0(6)*dt,zeta)
            call exact_ou_velocity(vpar,current_bg(7)**2,current_bg(6)*dt,zeta)
            call advance(-dt/2)
            phase_c=phase_c-dt*(bg0(4)+bg0(3)*(vc_old+vc)/2)
            if(anomalous_diffusion_coefficient>0) call diffuse(dt)
            ! Assign the exact endpoint when the clipped final step is taken.
            if(prefix_lag_nu>0.and.dt==target_time-t) then
                t=target_time
            else
                t=t+dt
            end if
            call sources()
            real_integral=real_integral+dt/2*(old_actual+actual_source)
            ref_integral=ref_integral+dt/2*(old_ref+ref_source)
            old_actual=actual_source;old_ref=ref_source
            if(.not.prefix_saved.and.t==prefix_time) then
                call ion_mu_deposit(bg0,wave,lambda,mu_ref,field0,v0,clight,real_integral,real_deposit)
                call ion_mu_deposit(bg0,wave,lambda,mu_ref,bg0(13),v0,clight,ref_integral,ref_deposit)
                prefix_result=real_deposit-ref_deposit
                prefix_saved=.true.
                prefix_checkpoint_count=prefix_checkpoint_count+1
            end if
        end do
        call ion_mu_deposit(bg0,wave,lambda,mu_ref,field0,v0,clight,real_integral,real_deposit)
        call ion_mu_deposit(bg0,wave,lambda,mu_ref,bg0(13),v0,clight,ref_integral,ref_deposit)
        result=real_deposit-ref_deposit
    contains
        subroutine diffuse(lag_step)
            ! Diffusion acts over positive lag; backward characteristics reverse
            ! only Hamiltonian advection. The real path uses the shared mesh
            ! displacement with reflection; the virtual flat cylinder receives
            ! the same draws in its fixed tangent frame (radius and helical phase).
            real(dp),intent(in) :: lag_step
            real(dp) :: delta_t,noise(3),increment(3)
            integer :: kicks,kick
            logical :: reflected
            kicks=max(1,ceiling(2*anomalous_diffusion_coefficient*lag_step/diffusion_step_cm**2))
            delta_t=lag_step/kicks
            do kick=1,kicks
                if(cylinder_only) then
                    call getran(0,noise(1));call getran(0,noise(2));noise(3)=0
                    ref_dr=ref_dr+sqrt(2*anomalous_diffusion_coefficient*delta_t)*noise(1)
                    ref_dphase=ref_dphase+bg0(2)*sqrt(2*anomalous_diffusion_coefficient*delta_t)*noise(2)
                else
                    if(debugging) call debug_state('kick_in',t,x,vpar,cell,face)
                    call anomalous_transport_displacement(x,cell,face,delta_t,vpar,vperp, &
                        anomalous_diffusion_coefficient,rho_bounds=diffusion_bounds,recover_lost=.false., &
                        reflected=reflected,random_vector=noise,time_sign=-1.0_dp)
                    if(debugging) call debug_state('kick_out',t,x,vpar,cell,face)
                    if(cell==-1) error stop 'reflected ion diffusion displacement lost'
                    if(.not.all(ieee_is_finite(x))) error stop 'invalid ion diffusion displacement'
                    call relocate_after_kick(x,vpar,vperp,mu,initialized,cell,face)
                    if(debugging) call debug_state('relocated',t,x,vpar,cell,face)
                    if(reflected) diffusion_reflections=diffusion_reflections+1
                    ref_dr=ref_dr+sqrt(2*delta_t)*dot_product(ref_radial,noise)
                    ref_dphase=ref_dphase+sqrt(2*delta_t)*dot_product(ref_helical,noise)
                end if
                diffusion_kicks=diffusion_kicks+1
            end do
            if(.not.cylinder_only) then
                rnow=flux_radius(x,cell);call sample_background(rnow,current_bg)
                max_excursion=max(max_excursion,abs(rnow-r0))
            end if
        end subroutine
        subroutine set_reference_frame(theta,Rlocal)
            real(dp),intent(in) :: theta,Rlocal
            call noise_frame(theta,Rlocal,ref_radial,ref_helical)
        end subroutine
        subroutine advance(step)
            real(dp),intent(in) :: step
            real(dp) :: saved_x(3),saved_vpar,saved_vperp,remaining
            integer :: saved_cell,saved_face,failure_unit
            if(cylinder_only) return
            saved_x=x;saved_vpar=vpar;saved_vperp=vperp;saved_cell=cell;saved_face=face
            if(debugging) call debug_state('step_in',t,x,vpar,cell,face)
            call orbit_timestep_gorilla(x,vpar,vperp,step,initialized,cell,face,remaining)
            if(debugging) call debug_state('step_out',t,x,vpar,cell,face)
            if(cell==-1) then
                open(newunit=failure_unit,file='ion_orbit_failure.dat',status='new')
                write(failure_unit,'(a)') 'GK_ION_ORBIT_FAILURE_V1'
                write(failure_unit,*) imu,radial_index,marker_index,saved_cell,saved_face
                write(failure_unit,'(7es26.17e3)') lambda,mu_ref,r0,t,step,remaining,vc
                write(failure_unit,'(5es26.17e3)') saved_x,saved_vpar,saved_vperp
                write(failure_unit,'(5es26.17e3)') x,vpar,vperp
                write(failure_unit,'(14es26.17e3)') current_bg
                close(failure_unit)
                error stop 'ion orbit lost; state captured; no respawn boundary selected'
            end if
            rnow=flux_radius(x,cell);call sample_background(rnow,current_bg)
            field=bmod_func(x-tetra_physics(cell)%x1,cell)
            max_mu_error=max(max_mu_error,abs(vperp*vperp/(2*field*mu)-1))
            max_excursion=max(max_excursion,abs(rnow-r0))
        end subroutine
        subroutine sources()
            complex(dp) :: phase(d)
            phase=exp(ii*(wave*(r0+ref_dr)+phase_c+ref_dphase))
            call ion_mu_source(bg0,wave,lambda,mu_ref,bg0(13),vc,phase,clight,ref_source)
            if(cylinder_only) then
                actual_source=ref_source
                ! Check mode deposits the diffused virtual path without subtraction.
                if(diffusion_reference_check) ref_source=0
                return
            end if
            alpha=mmode*theta_sfl(x)+nmode*n_field_periods_manual*x(2)
            field=bmod_func(x-tetra_physics(cell)%x1,cell)
            phase=exp(ii*(wave*rnow+alpha-alpha0))
            call ion_mu_source(current_bg,wave,lambda,mu_ref,field,vpar,phase,clight,actual_source)
        end subroutine
    end subroutine trace_marker

    subroutine trace_marker_paired(radial_index,marker_index,result,prefix_result,result0,prefix_result0)
        ! One marker, two real paths from the same injection: A with diffusion,
        ! B without. Both use the same collision draws on a common time grid;
        ! diffusion kicks come from a separate deterministic stream so that the
        ! collision stream stays synchronised. A's virtual cylinder receives the
        ! same kicks in its tangent frame; B's virtual cylinder is undiffused.
        use, intrinsic :: iso_fortran_env, only: int64
        integer,intent(in) :: radial_index,marker_index
        complex(dp),intent(out) :: result(d,d,9),prefix_result(d,d,9),result0(d,d,9),prefix_result0(d,d,9)
        real(dp) :: xa(3),xb(3),va,vb,vpa,vpb,v0,vc,vc_old,t,dt,rand2(2),theta0,alpha0,phase_c
        real(dp) :: bga(14),bgb(14),ra,rb,geometric,field0,mu,phase_rate,end_time,zeta,prefix_time,target_time
        real(dp) :: ref_dr,ref_dphase,ref_radial(3),ref_helical(3)
        complex(dp) :: srcA(d,3),srcB(d,3),refA(d,3),ref0(d,3),oA(d,3),oB(d,3),oRA(d,3),oR0(d,3)
        complex(dp) :: intA(d,3),intB(d,3),intRA(d,3),intR0(d,3),depA(d,d,9),depR(d,d,9)
        logical :: inita,initb,prefix_saved
        integer :: ca,fa,cb,fb,iteration
        integer(int64) :: kick_state
        call random_number(rand2)
        v0=bg0(7)*sqrt(-2*log(max(rand2(1),tiny(1.0_dp))))*cos(2*pi*rand2(2))
        va=v0;vb=v0;vc=v0;t=0;phase_c=0;mu=lambda*mu_ref;ref_dr=0;ref_dphase=0
        ra=r0;rb=r0;bga=bg0;bgb=bg0;field0=bg0(13)
        kick_state=int(seed,int64)*1000003_int64+int(imu,int64)*7919_int64 &
            +int(radial_index,int64)*104729_int64+int(marker_index,int64)*1299709_int64
        if(.not.cylinder_only) then
            call random_number(rand2);theta0=2*pi*rand2(1);geometric=r0
            vpa=sqrt(2*mu*bg0(13))
            do iteration=1,5
                xa=[R0_analytic_circ+geometric*cos(theta0),2*pi*rand2(2)/n_field_periods_manual,geometric*sin(theta0)]
                inita=.false.;ca=0;fa=0
                call orbit_timestep_gorilla(xa,va,vpa,0.0_dp,inita,ca,fa)
                if(ca==-1) error stop 'ion injection outside mesh'
                ra=flux_radius(xa,ca)
                if(abs(ra-r0)<1e-9_dp) exit
                geometric=geometric+r0-ra
            end do
            if(abs(ra-r0)>1e-8_dp) error stop 'ion injection flux label did not converge'
            field0=bmod_func(xa-tetra_physics(ca)%x1,ca)
            vpa=sqrt(2*mu*field0);inita=.false.
            call orbit_timestep_gorilla(xa,va,vpa,0.0_dp,inita,ca,fa)
            xb=xa;vpb=vpa;cb=ca;fb=fa;initb=.false.
            call orbit_timestep_gorilla(xb,vb,vpb,0.0_dp,initb,cb,fb)
            rb=ra;alpha0=mmode*theta_sfl(xa)+nmode*n_field_periods_manual*xa(2)
            call frame(theta0,xa(1))
        end if
        intA=0;intB=0;intRA=0;intR0=0
        call sources();oA=srcA;oB=srcB;oRA=refA;oR0=ref0
        end_time=lag_nu/bg0(6);prefix_time=prefix_lag_nu/bg0(6);prefix_saved=.false.
        prefix_result=0;prefix_result0=0
        do while(t<end_time)
            phase_rate=max(abs(bg0(4))+abs(bg0(3))*max(abs(vc),bg0(7)), &
                abs(bga(4))+abs(bga(3))*max(abs(va),bga(7)),abs(bgb(4))+abs(bgb(3))*max(abs(vb),bgb(7)))
            target_time=end_time
            if(.not.prefix_saved) target_time=prefix_time
            dt=min(nu_dt/max(bg0(6),bga(6),bgb(6)),phase_step/max(phase_rate,1.0_dp),target_time-t)
            vc_old=vc
            call step_both(-dt/2)
            call random_number(rand2)
            zeta=sqrt(-2*log(max(rand2(1),tiny(1.0_dp))))*cos(2*pi*rand2(2))
            call exact_ou_velocity(vc,bg0(7)**2,bg0(6)*dt,zeta)
            call exact_ou_velocity(va,bga(7)**2,bga(6)*dt,zeta)
            call exact_ou_velocity(vb,bgb(7)**2,bgb(6)*dt,zeta)
            call step_both(-dt/2)
            phase_c=phase_c-dt*(bg0(4)+bg0(3)*(vc_old+vc)/2)
            call diffuse_a(dt)
            if(dt==target_time-t) then
                t=target_time
            else
                t=t+dt
            end if
            call sources()
            intA=intA+dt/2*(oA+srcA);intB=intB+dt/2*(oB+srcB)
            intRA=intRA+dt/2*(oRA+refA);intR0=intR0+dt/2*(oR0+ref0)
            oA=srcA;oB=srcB;oRA=refA;oR0=ref0
            if(.not.prefix_saved.and.t==prefix_time) then
                call deposit(prefix_result,prefix_result0)
                prefix_saved=.true.;prefix_checkpoint_count=prefix_checkpoint_count+1
            end if
        end do
        call deposit(result,result0)
    contains
        subroutine deposit(withD,without)
            complex(dp),intent(out) :: withD(d,d,9),without(d,d,9)
            call ion_mu_deposit(bg0,wave,lambda,mu_ref,field0,v0,clight,intA,depA)
            call ion_mu_deposit(bg0,wave,lambda,mu_ref,bg0(13),v0,clight,intRA,depR)
            withD=depA-depR
            call ion_mu_deposit(bg0,wave,lambda,mu_ref,field0,v0,clight,intB,depA)
            call ion_mu_deposit(bg0,wave,lambda,mu_ref,bg0(13),v0,clight,intR0,depR)
            without=depA-depR
        end subroutine
        subroutine step_both(step)
            real(dp),intent(in) :: step
            if(cylinder_only) return
            call push(xa,va,vpa,inita,ca,fa,ra,bga,step)
            call push(xb,vb,vpb,initb,cb,fb,rb,bgb,step)
        end subroutine
        subroutine push(x,vpar,vperp,initialized,cell,face,rnow,current,step)
            real(dp),intent(inout) :: x(3),vpar,vperp,rnow,current(14)
            logical,intent(inout) :: initialized
            integer,intent(inout) :: cell,face
            real(dp),intent(in) :: step
            real(dp) :: saved_x(3),saved_vpar,saved_vperp,remaining
            integer :: saved_cell,saved_face,failure_unit
            saved_x=x;saved_vpar=vpar;saved_vperp=vperp;saved_cell=cell;saved_face=face
            call orbit_timestep_gorilla(x,vpar,vperp,step,initialized,cell,face,remaining)
            if(cell==-1) then
                open(newunit=failure_unit,file='ion_orbit_failure.dat',status='new')
                write(failure_unit,'(a)') 'GK_ION_ORBIT_FAILURE_V1'
                write(failure_unit,*) imu,radial_index,marker_index,saved_cell,saved_face
                write(failure_unit,'(7es26.17e3)') lambda,mu_ref,r0,t,step,remaining,vc
                write(failure_unit,'(5es26.17e3)') saved_x,saved_vpar,saved_vperp
                write(failure_unit,'(5es26.17e3)') x,vpar,vperp
                write(failure_unit,'(14es26.17e3)') current
                close(failure_unit)
                error stop 'ion orbit lost; state captured; no respawn boundary selected'
            end if
            rnow=flux_radius(x,cell);call sample_background(rnow,current)
            max_mu_error=max(max_mu_error,abs(vperp*vperp/(2*bmod_func(x-tetra_physics(cell)%x1,cell)*mu)-1))
            max_excursion=max(max_excursion,abs(rnow-r0))
        end subroutine
        subroutine diffuse_a(lag_step)
            real(dp),intent(in) :: lag_step
            real(dp) :: delta_t,noise(3),increment(3)
            integer :: kicks,kick,j
            logical :: reflected
            kicks=max(1,ceiling(2*anomalous_diffusion_coefficient*lag_step/diffusion_step_cm**2))
            delta_t=lag_step/kicks
            do kick=1,kicks
                do j=1,3
                    noise(j)=3.464102_dp*(uniform()-.5_dp)
                end do
                if(cylinder_only) then
                    ref_dr=ref_dr+sqrt(2*anomalous_diffusion_coefficient*delta_t)*noise(1)
                    ref_dphase=ref_dphase+bg0(2)*sqrt(2*anomalous_diffusion_coefficient*delta_t)*noise(2)
                else
                    call anomalous_transport_displacement(xa,ca,fa,delta_t,va,vpa, &
                        anomalous_diffusion_coefficient,rho_bounds=diffusion_bounds,recover_lost=.false., &
                        reflected=reflected,noise_in=noise,time_sign=-1.0_dp)
                    if(ca==-1) error stop 'reflected ion diffusion displacement lost'
                    if(.not.all(ieee_is_finite(xa))) error stop 'invalid ion diffusion displacement'
                    call relocate_after_kick(xa,va,vpa,mu,inita,ca,fa)
                    if(reflected) diffusion_reflections=diffusion_reflections+1
                    ref_dr=ref_dr+sqrt(2*delta_t)*dot_product(ref_radial,noise)
                    ref_dphase=ref_dphase+sqrt(2*delta_t)*dot_product(ref_helical,noise)
                end if
                diffusion_kicks=diffusion_kicks+1
            end do
            if(.not.cylinder_only) then
                ra=flux_radius(xa,ca);call sample_background(ra,bga)
                max_excursion=max(max_excursion,abs(ra-r0))
            end if
        end subroutine
        real(dp) function uniform()
            ! splitmix64; independent of the intrinsic collision stream.
            integer(int64) :: z
            kick_state=kick_state+(-7046029254386353131_int64)
            z=kick_state
            z=ieor(z,ishft(z,-30))*(-4658895280553007687_int64)
            z=ieor(z,ishft(z,-27))*(-7723592293110705685_int64)
            z=ieor(z,ishft(z,-31))
            uniform=(real(ishft(z,-11),dp)+.5_dp)*2.0_dp**(-53)
        end function
        subroutine frame(theta,Rlocal)
            real(dp),intent(in) :: theta,Rlocal
            call noise_frame(theta,Rlocal,ref_radial,ref_helical)
        end subroutine
        subroutine sources()
            complex(dp) :: phase(d)
            real(dp) :: alpha
            phase=exp(ii*(wave*(r0+ref_dr)+phase_c+ref_dphase))
            call ion_mu_source(bg0,wave,lambda,mu_ref,bg0(13),vc,phase,clight,refA)
            phase=exp(ii*(wave*r0+phase_c))
            call ion_mu_source(bg0,wave,lambda,mu_ref,bg0(13),vc,phase,clight,ref0)
            if(cylinder_only) then
                srcA=refA;srcB=ref0;return
            end if
            alpha=mmode*theta_sfl(xa)+nmode*n_field_periods_manual*xa(2)
            phase=exp(ii*(wave*ra+alpha-alpha0))
            call ion_mu_source(bga,wave,lambda,mu_ref,bmod_func(xa-tetra_physics(ca)%x1,ca),va,phase,clight,srcA)
            alpha=mmode*theta_sfl(xb)+nmode*n_field_periods_manual*xb(2)
            phase=exp(ii*(wave*rb+alpha-alpha0))
            call ion_mu_source(bgb,wave,lambda,mu_ref,bmod_func(xb-tetra_physics(cb)%x1,cb),vb,phase,clight,srcB)
        end subroutine
    end subroutine trace_marker_paired

    subroutine debug_state(tag,t,x,vpar,cell,face)
        ! Signed (unnormalised) face distances and the pusher's roundoff bound.
        character(*),intent(in) :: tag
        real(dp),intent(in) :: t,x(3),vpar
        integer,intent(in) :: cell,face
        real(dp) :: dist(4),bound(4),scale(3)
        integer :: l
        dist=0;bound=0
        if(cell>0) then
            scale=2*abs(x)+abs(tetra_physics(cell)%x1)
            do l=1,4
                dist(l)=sum(tetra_physics(cell)%anorm(:,l)*(x-tetra_physics(cell)%x1))
                bound(l)=sum(abs(tetra_physics(cell)%anorm(:,l))*scale)+abs(tetra_physics(cell)%dist_ref_vec(l))
                if(l==1) dist(l)=dist(l)+tetra_physics(cell)%dist_ref
                if(l==1) bound(l)=bound(l)+abs(tetra_physics(cell)%dist_ref)
            end do
            bound=8*epsilon(1.0_dp)*bound
        end if
        write(debug_unit,'(a10,es24.15,3es24.15,es24.15,2i10,8es12.3)') tag,t,x,vpar,cell,face,dist,bound
        flush(debug_unit)
    end subroutine

    subroutine relocate_after_kick(x,vpar,vperp,mu,initialized,cell,face)
        ! The straight-line diffusion displacement can end a numerical hair
        ! outside the face it reports as the entry face. The pusher never exits
        ! through its entry face, so a marker there whose Hamiltonian velocity
        ! points outward cannot proceed. Geometrically it lies in the neighbour
        ! across that face: hand it over with the pusher's own mesh logic
        ! (including periodic shifts). Position, time, velocity and mu unchanged.
        ! (A velocity-disambiguated cell search does not suffice: it can return
        ! the same cell when E x B dominates the face-normal velocity.)
        real(dp),intent(inout) :: x(3),vpar,vperp
        real(dp),intent(in) :: mu
        logical,intent(inout) :: initialized
        integer,intent(inout) :: cell,face
        real(dp) :: distance
        integer :: neighbour,iper
        if(face==0) return
        distance=sum(tetra_physics(cell)%anorm(:,face)*(x-tetra_physics(cell)%x1))
        if(face==1) distance=distance+tetra_physics(cell)%dist_ref
        ! Within the kick's face-convergence tolerance (1e-10 dist_ref), the
        ! displacement routine already chose the cell from the traced-time normal
        ! velocity; a handover there would undo that choice. Only a marker
        ! outside its entry face beyond that tolerance belongs to the neighbour.
        if(distance>=-1.0e-10_dp*abs(tetra_physics(cell)%dist_ref)) return
        call pusher_handover2neighbour(cell,neighbour,face,x,iper)
        if(neighbour<1) error stop 'ion diffusion handover across the mesh boundary'
        cell=neighbour
        diffusion_relocations=diffusion_relocations+1
    end subroutine

    subroutine noise_frame(theta,Rlocal,radial,helical)
        ! Noise-space vectors of the flat virtual cylinder at injection: the
        ! radial displacement is radial.xi*sqrt(2 dt) and the helical phase
        ! helical.xi*sqrt(2 dt), for the same three unit-variance draws xi as the
        ! real kick. Start from the injection Cholesky factor (maximal correlation
        ! with the real path), then orthonormalise in noise space so the virtual
        ! covariance is exactly D and D*ks^2 with zero cross term. The (R,phi,Z)
        ! factor carries roundoff from the 1/R^2 metric range; the check only
        ! rejects a gross frame error.
        real(dp),intent(in) :: theta,Rlocal
        real(dp),intent(out) :: radial(3),helical(3)
        real(dp) :: a,c,htheta,hphi,norm,h(3),s,co,factor(3,3),rr,hh
        a=mmode/r0;c=nmode*n_field_periods_manual/R0_analytic_circ
        htheta=(a*bg0(3)-c*bg0(2))/(a*a+c*c);hphi=(c*bg0(3)+a*bg0(2))/(a*a+c*c)
        norm=sqrt(htheta*htheta+hphi*hphi);htheta=htheta/norm;hphi=hphi/norm
        s=sin(theta);co=cos(theta);h=[-htheta*s,hphi/Rlocal,htheta*co]
        call compute_diffusion_cholesky(h,Rlocal,anomalous_diffusion_coefficient,factor)
        radial=matmul(transpose(factor),[co,0.0_dp,s])
        helical=matmul(transpose(factor),bg0(2)*[-hphi*s,-htheta*Rlocal,hphi*co])
        rr=sum(radial**2);hh=sum(helical**2)
        if(abs(rr/anomalous_diffusion_coefficient-1)>1e-3_dp.or. &
            abs(hh/(anomalous_diffusion_coefficient*bg0(2)**2)-1)>1e-3_dp.or. &
            abs(dot_product(radial,helical))>1e-3_dp*sqrt(rr*hh)) &
            error stop 'ion diffusion reference frame is not the perpendicular tangent frame'
        radial=radial*sqrt(anomalous_diffusion_coefficient/rr)
        helical=helical-dot_product(helical,radial)/anomalous_diffusion_coefficient*radial
        helical=helical*abs(bg0(2))*sqrt(anomalous_diffusion_coefficient/sum(helical**2))
    end subroutine

    subroutine sample_background(radius,result)
        real(dp),intent(in) :: radius
        real(dp),intent(out) :: result(14)
        real(dp) :: index,f
        integer :: lo,hi
        index=modulo(radius-bg(1,1),L)*N/L
        lo=floor(index)+1;hi=mod(lo,N)+1;f=index-floor(index)
        result=(1-f)*bg(lo,:)+f*bg(hi,:);result(1)=radius
    end subroutine
    real(dp) function flux_radius(position,cell) result(radius)
        real(dp),intent(in) :: position(3)
        integer,intent(in) :: cell
        real(dp) :: psi,f
        integer :: lo,hi,mid
        psi=abs(tetra_physics(cell)%Aphi1+sum(tetra_physics(cell)%gAphi*(position-tetra_physics(cell)%x1)))
        if(psi<flux_psi(1).or.psi>flux_psi(flux_points)) error stop 'ion flux radius outside mapping'
        lo=1;hi=flux_points
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
        real(dp),intent(in) :: position(3)
        real(dp) :: angle,eps,radius
        radius=sqrt((position(1)-R0_analytic_circ)**2+position(3)**2)
        angle=atan2(position(3),position(1)-R0_analytic_circ);eps=radius/R0_analytic_circ
        theta=2*atan2(sqrt(1-eps)*sin(angle/2),sqrt(1+eps)*cos(angle/2))
    end function
    subroutine write_response(path,result)
        character(*),intent(in) :: path
        complex(dp),intent(in) :: result(:,:,:)
        integer :: u,row,col,k,j
        if(.not.all(ieee_is_finite(real(result))).or..not.all(ieee_is_finite(aimag(result)))) &
            error stop 'nonfinite ion orbit response'
        open(newunit=u,file=path,status='new')
        write(u,'(a)') 'GK_ION_RESPONSE_V1'
        write(u,*) M,N,model,mmode,nmode,species
        write(u,'(3es26.17e3)') L,rm,clight
        do j=1,N;write(u,'(14es26.17e3)') bg(j,:);end do
        do k=1,4
            do col=1,d
                do row=1,d
                    write(u,'(2es26.17e3)') real(result(row,col,k)),aimag(result(row,col,k))
                end do
            end do
        end do
        close(u)
    end subroutine
end module

program gorilla_orbit_ion_response
    use orbit_ion_response_driver_mod, only: run_orbit_ion_response
    implicit none
    call run_orbit_ion_response()
end program
