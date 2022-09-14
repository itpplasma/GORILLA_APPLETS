!
    module gorilla_applets_settings_mod
!
        implicit none
!
        private
!
        !Fluxtube volume OR mono-energetic radial transport coefficient OR numerical diffusion coefficient
        integer, public, protected :: i_option
!
        !Options for pre-computation of fluxtube volume
        character(1024),public,protected :: filename_fluxtv_precomp
        double precision, public, protected :: start_pos_x1,start_pos_x2,start_pos_x3,t_step_fluxtv,energy_eV_fluxtv
        integer(kind=8), public, protected :: nt_steps_fluxtv
        character(1024),public,protected :: filename_fluxtv_load
!
        public :: load_gorilla_applets_inp

        NAMELIST /GORILLA_APPLETS_NML/ i_option, filename_fluxtv_precomp, start_pos_x1, start_pos_x2, start_pos_x3, t_step_fluxtv, &
            & energy_eV_fluxtv,nt_steps_fluxtv,filename_fluxtv_load
!
    contains
!            
        subroutine load_gorilla_applets_inp()
!
            open(unit=11, file='gorilla_applets.inp', status='unknown')
            read(11,nml=GORILLA_APPLETS_NML)
            close(11)

            print *,'GORILLA_APPLETS: Loaded input data from input file'
!            
        end subroutine load_gorilla_applets_inp

    end module gorilla_applets_settings_mod
!    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
