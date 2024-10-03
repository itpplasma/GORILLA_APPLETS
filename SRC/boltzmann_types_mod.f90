module boltzmann_types_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    type moment_specs_t
    logical :: boole_squared_moments
    integer, dimension(4) :: moments_selector
    integer :: n_moments
    integer :: n_triangles
    integer :: n_fourier_modes
    end type moment_specs_t

    type output_t
    real(dp), dimension(:), allocatable :: prism_volumes
    real(dp), dimension(:), allocatable :: refined_prism_volumes
    real(dp), dimension(:), allocatable :: electric_potential
    real(dp), dimension(:), allocatable :: boltzmann_density
    complex(dp), dimension(:,:), allocatable :: tetr_moments
    complex(dp), dimension(:,:), allocatable :: prism_moments
    complex(dp), dimension(:,:), allocatable :: prism_moments_squared
    complex(dp), dimension(:,:,:), allocatable :: moments_in_frequency_space
    end type output_t

    type boole_writing_data_t
    logical :: vertex_indices = .false.
    logical :: vertex_coordinates = .false.
    logical :: prism_volumes = .false.
    logical :: refined_prism_volumes = .true.
    logical :: boltzmann_density = .true.
    logical :: electric_potential = .false.
    logical :: moments = .false.
    logical :: fourier_moments = .true.
    end type boole_writing_data_t

    type filenames_t
    character(len=100) :: exit_times
    character(len=100) :: remaining_particles
    character(len=100) :: poincare_maps
    character(len=100) :: prism_moments
    character(len=100) :: prism_moments_summed_squares
    character(len=100) :: vertex_coordinates
    character(len=100) :: vertex_indices
    character(len=100) :: prism_volumes
    character(len=100) :: fourier_moments
    character(len=100) :: refined_prism_volumes
    character(len=100) :: electric_potential
    character(len=100) :: boltzmann_density
    character(len=100) :: divertor_intersections
    character(len=100) :: tetr_moments
    end type filenames_t

    type counter_t
    integer :: lost_particles = 0
    integer :: lost_inside = 0
    integer :: lost_outside = 0
    integer :: tetr_pushings = 0
    integer :: phi_0_mappings = 0
    integer :: integration_steps = 0
    end type counter_t

    type poloidal_flux_t
    real(dp), dimension(:), allocatable :: particle
    real(dp) :: min
    real(dp) :: max
    end type poloidal_flux_t

    type collisions_t
    real(dp), dimension(:,:), allocatable :: dens_mat
    real(dp), dimension(:,:), allocatable :: temp_mat
    real(dp), dimension(:,:), allocatable :: vpar_mat
    real(dp), dimension(:,:), allocatable :: efcolf_mat
    real(dp), dimension(:,:), allocatable :: velrat_mat
    real(dp), dimension(:,:), allocatable :: enrat_mat
    real(dp), dimension(:), allocatable :: mass_num
    real(dp), dimension(:), allocatable :: charge_num
    real(dp), dimension(:), allocatable :: dens
    real(dp), dimension(:), allocatable :: temp
    end type collisions_t

    type boltzmann_input_t
    real(dp) :: time_step
    real(dp) :: energy_eV
    real(dp) :: n_particles
    real(dp) :: density
    logical  :: boole_squared_moments
    logical  :: boole_point_source
    logical  :: boole_collisions
    logical  :: boole_precalc_collisions
    logical  :: boole_refined_sqrt_g
    logical  :: boole_boltzmann_energies
    logical  :: boole_linear_density_simulation
    logical  :: boole_antithetic_variate
    logical  :: boole_linear_temperature_simulation
    integer  :: i_integrator_type
    integer  :: seed_option
    end type boltzmann_input_t

end module boltzmann_types_mod