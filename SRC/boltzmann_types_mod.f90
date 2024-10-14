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

    type(moment_specs_t) :: moment_specs

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

    type(output_t) :: output

    type start_t
    real(dp), dimension(:,:), allocatable :: x
    real(dp), dimension(:), allocatable :: pitch
    real(dp), dimension(:), allocatable :: energy
    real(dp), dimension(:), allocatable :: weight
    real(dp), dimension(:), allocatable :: jperp
    end type start_t

    type(start_t) :: start

    type counter_t
    integer :: lost_particles = 0
    integer :: lost_inside = 0
    integer :: lost_outside = 0
    integer :: tetr_pushings = 0
    integer :: phi_0_mappings = 0
    integer :: integration_steps = 0
    end type counter_t

    type(counter_t) :: counter

    type poloidal_flux_t
    real(dp) :: min
    real(dp) :: max
    end type poloidal_flux_t

    type(poloidal_flux_t) :: pflux

    type collisions_t
    real(dp), dimension(:,:,:), allocatable :: randcol
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
    integer :: randcoli = int(1.0d5)
    integer :: n !number of background species
    real(dp) :: maxcol = 0
    end type collisions_t

    type(collisions_t) :: c

    type input_t
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
    integer  :: num_particles
    logical  :: boole_write_vertex_indices
    logical  :: boole_write_vertex_coordinates
    logical  :: boole_write_prism_volumes
    logical  :: boole_write_refined_prism_volumes
    logical  :: boole_write_boltzmann_density
    logical  :: boole_write_electric_potential
    logical  :: boole_write_moments
    logical  :: boole_write_fourier_moments
    logical  :: boole_write_exit_data
    logical  :: boole_divertor_intersections !Used in divertor_heat_loads
    logical  :: boole_poincare_mappings !Used in divertor_heat_loads
    integer  :: num_poincare_mappings !Used in divertor_heat_loads
    end type input_t

    type(input_t) :: u

    type filenames_t
    character(len=100) :: exit_times
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
    character(len=100) :: exit_data
    end type filenames_t

    type(filenames_t) :: filenames

    type exit_data_t
    integer, dimension(:), allocatable :: lost
    real(dp), dimension(:), allocatable :: t_confined
    real(dp), dimension(:,:), allocatable :: x
    real(dp), dimension(:), allocatable :: vpar
    real(dp), dimension(:), allocatable :: vperp
    integer, dimension(:), allocatable :: integration_step
    integer, dimension(:), allocatable :: phi_0_mappings
    end type exit_data_t

    type(exit_data_t) :: exit_data

    type grid_t
    real(dp) :: amin
    real(dp) :: amax
    real(dp) :: cmin
    real(dp) :: cmax
    integer  :: ind_a
    integer  :: ind_b
    integer  :: ind_c
    end type grid_t

    type(grid_t) :: g

    type boole_t
    logical :: initialized
    logical :: lost
    logical :: exit !used in divertor heat loads mod
    end type boole_t

    type time_t
    real(dp) :: step
    real(dp) :: remain
    real(dp) :: confined
    end type time_t

end module boltzmann_types_mod