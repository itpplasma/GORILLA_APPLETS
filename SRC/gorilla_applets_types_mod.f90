module gorilla_applets_types_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use marker_distribution_mod, only: distribution_3d_t, distribution_1d_t

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
    real(dp), dimension(:), allocatable :: radial_flux
    real(dp), dimension(:), allocatable :: weighted_radial_flux
    complex(dp), dimension(:,:,:), allocatable :: tetr_moments
    complex(dp), dimension(:,:,:), allocatable :: prism_moments
    complex(dp), dimension(:,:,:), allocatable :: prism_moments_squared
    complex(dp), dimension(:,:,:), allocatable :: moments_in_frequency_space
    !In moments_in_frequency_space, third dimension is not particle species but fourier mode
    end type output_t

    type(output_t) :: output

    type start_t
    real(dp), dimension(:,:,:), allocatable :: x
    real(dp), dimension(:,:),   allocatable :: pitch
    real(dp), dimension(:,:),   allocatable :: energy
    real(dp), dimension(:,:),   allocatable :: jperp
    real(dp), dimension(:),     allocatable :: particle_charge !used in self consistent electric field computation
    real(dp), dimension(:),     allocatable :: particle_mass !used in self consistent electric field computation
    real(dp), dimension(:),     allocatable :: cm_over_e !used in self consistent electric field computation
    real(dp), dimension(:),     allocatable :: v0 !This is the thermal velocity, not the particle velocity;
                                                  !used in self consistent electric field computation
    real(dp), dimension(:),     allocatable :: t !total tracing time
    real(dp)                                :: epsilon_max
    logical, dimension(:,:),    allocatable :: lost
    type(distribution_3d_t) :: dist_position  ! starting distribution in position space
    type(distribution_1d_t) :: dist_energy    ! starting distribution in energy
    type(distribution_1d_t) :: dist_lambda    ! starting distribution in pitch angle
    end type start_t

    type(start_t) :: start

    type counter_t
    integer :: lost_particles = 0
    integer :: lost_inside = 0
    integer :: tetr_pushings = 0
    integer :: phi_0_mappings = 0
    integer :: integration_steps = 0
    end type counter_t

    type(counter_t) :: counter

    type flux_t
    real(dp) :: poloidal_min
    real(dp) :: poloidal_max
    end type flux_t

    type(flux_t) :: flux

    type collisions_t
    real(dp), dimension(:,:,:,:), allocatable :: randcol
    real(dp), dimension(:,:), allocatable :: dens_mat
    real(dp), dimension(:,:), allocatable :: temp_mat
    real(dp), dimension(:,:), allocatable :: vpar_mat
    real(dp), dimension(:,:), allocatable :: efcolf_mat
    real(dp), dimension(:,:), allocatable :: velrat_mat
    real(dp), dimension(:,:), allocatable :: enrat_mat
    real(dp), dimension(:), allocatable :: mass
    real(dp), dimension(:), allocatable :: charge_num
    real(dp), dimension(:), allocatable :: dens
    real(dp), dimension(:), allocatable :: temp
    integer :: randcoli = int(1.0d5)
    integer :: n !number of background species
    real(dp) :: maxcol = 0
    real(dp) :: weight_factor
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
    logical  :: boole_monoenergetic
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
    logical  :: boole_write_grid_data
    logical  :: boole_preserve_energy_and_momentum_during_collisions = .false.
    integer  :: n_background_density_updates = 0
    logical  :: boole_static_ne !used in self consistent electric field computation
    integer  :: n_species = 1 !used in self consistent electric field computation
    integer  :: tracer_species = 1
    integer  :: n_electric_potential_updates !used in self consistent electric field computation
    integer  :: update_dimension = 1 !used in self consistent electric field computation
    logical  :: boole_calc_diffusion_coefficient !used in anomalous transport
    integer  :: i_scan_option = 0 !used in anomalous transport: 0=no scan, 1=scan over eps_Phi, 2=scan over n2, 3=scan over n3
    real(dp) :: anomalous_diffusion_coefficient   !anomalous diffusion coefficient in cm^2/s
    logical  :: boole_divertor_intersection !Used in divertor_heat_loads
    logical  :: boole_poincare_plot !Used in divertor_heat_loads
    logical  :: boole_eliminate_particles_outside_flux !Used in helical_core
    real(dp) :: flux_threshold_for_elimination !Used in helical_core (0=axis, 1=boundary)
    logical  :: boole_delta_f = .false. !Used in helical_core: use delta-f method for weights
    integer  :: n_poincare_mappings !Used in divertor_heat_loads
    integer  :: n_mappings_ignored !Used in divertor_heat_loads
    integer  :: collision_operator = 4
    logical  :: boole_custom_background = .false.
    real(dp), dimension(2) :: background_density_cm3 = 0.0_dp
    real(dp), dimension(2) :: background_temperature_eV = 0.0_dp
    logical  :: boole_custom_source_profile = .false.
    real(dp) :: density_log_gradient_per_s = 0.0_dp
    real(dp) :: density_profile_reference_s = 0.5_dp
    real(dp) :: source_profile_reference_s = 0.5_dp
    real(dp) :: source_profile_width = 0.1_dp
    real(dp) :: lambda !Used in divertor_heat_loads
    end type input_t

    type(input_t) :: in

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
    character(len=100) :: grid_data
    end type filenames_t

    type(filenames_t) :: filenames

    type exit_data_t
    integer, dimension(:,:), allocatable :: lost
    real(dp), dimension(:,:), allocatable :: t_confined
    real(dp), dimension(:,:,:), allocatable :: x
    real(dp), dimension(:,:), allocatable :: vpar
    real(dp), dimension(:,:), allocatable :: vperp
    integer, dimension(:,:), allocatable :: integration_step
    integer, dimension(:,:), allocatable :: phi_0_mappings
    real(dp), dimension(:,:), allocatable :: flux_surface
    end type exit_data_t

    type(exit_data_t) :: exit_data

    type grid_t
    real(dp) :: amin
    real(dp) :: amax
    real(dp) :: cmin
    real(dp) :: cmax
    real(dp) :: raxis
    real(dp) :: zaxis
    real(dp) :: dist_from_o_point_within_grid !minimum radius of circle centered at magnetic axis completely contained in the grid
    real(dp) :: total_volume
    integer  :: ind_a
    integer  :: ind_b
    integer  :: ind_c
    integer, dimension(:,:), allocatable :: vertices_per_flux_surface !used in self_consistent_electric_field_mod
    integer, dimension(:,:), allocatable :: prisms_per_flux_tube !used in self_consistent_electric_field_mod
    end type grid_t

    type(grid_t) :: g

    type electric_potential_t
    real(dp), dimension(:), allocatable :: rho_prism
    real(dp), dimension(:), allocatable :: rho_flux_layer
    real(dp), dimension(:), allocatable :: rho_vert
    real(dp), dimension(:), allocatable :: phi_elec_from_rho
    real(dp), dimension(:), allocatable :: average_abs_phi_elec_from_rho
    real(dp), dimension(:), allocatable :: total_tracing_time
    real(dp), dimension(:), allocatable :: s_shell_volumes
    real(dp) :: mean_abs_rho_at_first_update
    end type electric_potential_t

    type(electric_potential_t) :: ep

    type delta_s_delta_s_squared_t
    real(dp), dimension(:), allocatable    :: delta_s
    real(dp), dimension(:), allocatable    :: delta_s_squared
    logical, dimension(:), allocatable     :: boole_large_distance
    real(dp), dimension(:), allocatable    :: time
    real(dp), dimension(:,:), allocatable  :: f_v
    integer                                :: k, j
    real(dp)                               :: s0
    integer, dimension(:), allocatable     :: check
    integer                                :: n_particles
    real(dp)                               :: temperature
    real(dp)                               :: integral_for_weight_normalisation
    end type delta_s_delta_s_squared_t

    type(delta_s_delta_s_squared_t) :: s

    type diffusion_coefficients 
    real(dp), dimension(:), allocatable :: A
    real(dp), dimension(:), allocatable :: A_from_first_run
    real(dp), dimension(:), allocatable :: B
    real(dp), dimension(:), allocatable :: grad_A
    real(dp), dimension(:), allocatable :: grad_B
    real(dp), dimension(:), allocatable :: s_vertices
    real(dp), dimension(4)              :: polynomial_coefficients_for_B
    real(dp), dimension(3)              :: polynomial_coefficients_for_A
    end type diffusion_coefficients

    type(diffusion_coefficients) dc

    type one_d_t
    real(dp), dimension(:,:), allocatable :: densities
    logical                               :: boole_print_densities
    end type one_d_t

    type(one_d_t) one_d

    type particle_status_t
    logical :: initialized
    logical :: lost
    logical :: exit !used in divertor heat loads mod
    end type particle_status_t

    type time_t
    real(dp) :: step
    real(dp) :: remain
    real(dp) :: confined
    real(dp) :: step_anomalous_transport
    end type time_t

    type weights_t
    real(dp), dimension(:,:), allocatable :: w         ! current weight (particles x species)
    real(dp), dimension(:,:), allocatable :: original  ! original weight, used in helical_core for delta-f damping
    end type weights_t

    type(weights_t) :: weights

    real(dp) :: maximum_s

end module gorilla_applets_types_mod
