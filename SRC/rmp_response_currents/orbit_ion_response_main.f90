module orbit_ion_response_driver_mod
    ! Large-aspect-ratio ion response on actual equilibrium tetrahedral paths.
    ! Energy is closed at each GLOBAL magnetic moment before mu integration.
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use omp_lib, only: omp_get_wtime
    use cylindrical_ion_mu_response_mod, only: ion_mu_source,ion_mu_deposit,assemble_cylindrical_ion_mu
    use cylindrical_ion_response_mod, only: perpendicular_rule
    use cylindrical_electron_response_mod, only: close_energy_response
    use orbit_timestep_gorilla_mod, only: initialize_gorilla,orbit_timestep_gorilla
    use gorilla_settings_mod, only: load_gorilla_inp,ipusher,coord_system, &
        boole_helical_pert,boole_pert_from_mephit,boole_axi_noise_vector_pot, &
        boole_axi_noise_elec_pot,boole_non_axi_noise_vector_pot,boole_strong_electric_field,boole_magnetic_drift
    use tetra_grid_settings_mod, only: load_tetra_grid_inp,R0_analytic_circ,n_field_periods_manual,grid_kind
    use tetra_physics_mod, only: tetra_physics,particle_mass,particle_charge
    use supporting_functions_mod, only: bmod_func,energy_tot_func
    use field_analytic_circ_mod, only: psi_pol_analytic_circ
    use collis_ions, only: exact_ou_velocity
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
    namelist /ion_orbit_response/ radial_samples,markers,batches,n_mu,seed,lag_nu,nu_dt,phase_step, &
        mu_reference_scale,cylinder_only,prefix_lag_nu,radial_sample_shift_cm
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
        open(newunit=paired_unit,file='ion_paired_mu.dat',status='new')
        write(paired_unit,'(a)') '# node batch window channel column row real imag'
        write(paired_unit,'(a)') '# window 0=bare, 1=prefix correction, 2=total correction, 3=exact closed'
    end if
    open(newunit=progress,file='ion_orbit_progress.dat',status='new')
    write(progress,'(a)') '# mu_node radial_sample elapsed_seconds'
    open(newunit=closure_unit,file='ion_mu_closure.dat',status='new')
    write(closure_unit,'(a)') '# node lambda weight mean_closure_residual mean_rcond'
    start_time=omp_get_wtime()
    do imu=1,n_mu
        lambda=nodes(imu)
        call assemble_cylindrical_ion_mu(bg,M,L,clight,lambda,mu_ref,model,exact,bare,refinement)
        integrated_exact=integrated_exact+weights(imu)*exact
        if(model==1) call close_energy_response(bare,cylinder_closed,closure_error,rcond)
        corrections=0
        if(prefix_lag_nu>0) prefix_corrections=0
        do ir=1,radial_samples
            r0=bg(1,1)+(ir-.5_dp)*L/radial_samples+radial_sample_shift_cm
            call sample_background(r0,bg0)
            projection=exp(-ii*wave*r0)/real(radial_samples*markers/batches,dp)
            do im=1,markers
                call trace_marker(ir,im,sample_result,prefix_sample)
                ib=mod(im-1,batches)+1
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
            call write_mu_channels(0,0,bare)
            call write_mu_channels(0,3,exact)
            do ib=1,batches
                call write_mu_channels(ib,1,prefix_corrections(:,:,:,ib))
                call write_mu_channels(ib,2,corrections(:,:,:,ib))
                if(model==1) then
                    call close_energy_response(bare+prefix_corrections(:,:,:,ib),closed,closure_error,rcond)
                    prefix_batches(:,:,:,ib)=prefix_batches(:,:,:,ib)+weights(imu)*(exact+closed-cylinder_closed)
                    deallocate(closed)
                else
                    prefix_batches(:,:,:,ib)=prefix_batches(:,:,:,ib)+weights(imu)*(exact+prefix_corrections(:,:,1:4,ib))
                end if
            end do
            flush(paired_unit)
        end if
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
    if(prefix_lag_nu>0) then
        close(prefix_closure_unit);close(paired_unit)
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
    if(prefix_lag_nu>0) then
        if(prefix_checkpoint_count/=n_mu*radial_samples*markers) error stop 'missing marker prefix checkpoint'
        write(diagnostics,'(a,i12)') 'prefix_checkpoint_count ',prefix_checkpoint_count
    end if
    close(diagnostics)
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
        call random_number(rand2)
        v0=bg0(7)*sqrt(-2*log(max(rand2(1),tiny(1.0_dp))))*cos(2*pi*rand2(2))
        vpar=v0;vc=v0;t=0;phase_c=0;mu=lambda*mu_ref
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
        subroutine advance(step)
            real(dp),intent(in) :: step
            real(dp) :: saved_x(3),saved_vpar,saved_vperp,remaining
            integer :: saved_cell,saved_face,failure_unit
            if(cylinder_only) return
            saved_x=x;saved_vpar=vpar;saved_vperp=vperp;saved_cell=cell;saved_face=face
            call orbit_timestep_gorilla(x,vpar,vperp,step,initialized,cell,face,remaining)
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
            phase=exp(ii*(wave*r0+phase_c))
            call ion_mu_source(bg0,wave,lambda,mu_ref,bg0(13),vc,phase,clight,ref_source)
            if(cylinder_only) then
                actual_source=ref_source;return
            end if
            alpha=mmode*theta_sfl(x)+nmode*n_field_periods_manual*x(2)
            field=bmod_func(x-tetra_physics(cell)%x1,cell)
            phase=exp(ii*(wave*rnow+alpha-alpha0))
            call ion_mu_source(current_bg,wave,lambda,mu_ref,field,vpar,phase,clight,actual_source)
        end subroutine
    end subroutine trace_marker

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
