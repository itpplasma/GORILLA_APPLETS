module utils_anomalous_transport_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine read_anomalous_transport_inp_into_type

    use gorilla_applets_types_mod, only: in

    real(dp) :: time_step,energy_eV,n_particles, density
    logical :: boole_squared_moments, boole_point_source, boole_collisions, boole_precalc_collisions, boole_refined_sqrt_g, &
               boole_boltzmann_energies, boole_linear_density_simulation, boole_antithetic_variate, &
               boole_linear_temperature_simulation, boole_write_vertex_indices, boole_write_vertex_coordinates, &
               boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
               boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
               boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, boole_static_ne
    integer :: i_integrator_type, seed_option, n_electric_potential_updates, update_dimension, n_species

    integer :: s_inp_unit

    !Namelist for anomalous transport input (uses same namelist name for compatibility)
    NAMELIST /self_consistent_ef_nml/ time_step,energy_eV,n_particles,boole_squared_moments,boole_point_source,boole_collisions, &
    & boole_precalc_collisions,density,boole_refined_sqrt_g,boole_boltzmann_energies, boole_linear_density_simulation, &
    & boole_antithetic_variate,boole_linear_temperature_simulation,i_integrator_type,seed_option, boole_write_vertex_indices, &
    & boole_write_vertex_coordinates, boole_write_prism_volumes, boole_write_refined_prism_volumes, boole_write_boltzmann_density, &
    & boole_write_electric_potential, boole_write_moments, boole_write_fourier_moments, boole_write_exit_data, &
    & boole_write_grid_data, boole_preserve_energy_and_momentum_during_collisions, n_electric_potential_updates, update_dimension, &
    & n_species, boole_static_ne

    open(newunit = s_inp_unit, file='anomalous_transport.inp', status='unknown')
    read(s_inp_unit,nml=self_consistent_ef_nml)
    close(s_inp_unit)

    in%time_step = time_step
    in%energy_eV = energy_eV
    in%n_particles = n_particles
    in%density = density
    in%boole_squared_moments = boole_squared_moments
    in%boole_point_source = boole_point_source
    in%boole_collisions = boole_collisions
    in%boole_precalc_collisions = boole_precalc_collisions
    in%boole_refined_sqrt_g = boole_refined_sqrt_g
    in%boole_boltzmann_energies = boole_boltzmann_energies
    in%boole_linear_density_simulation = boole_linear_density_simulation
    in%boole_antithetic_variate = boole_antithetic_variate
    in%boole_linear_temperature_simulation = boole_linear_temperature_simulation
    in%i_integrator_type = i_integrator_type
    in%seed_option = seed_option
    in%num_particles = int(n_particles)
    in%boole_write_vertex_indices = boole_write_vertex_indices
    in%boole_write_vertex_coordinates = boole_write_vertex_coordinates
    in%boole_write_prism_volumes = boole_write_prism_volumes
    in%boole_write_refined_prism_volumes = boole_write_refined_prism_volumes
    in%boole_write_boltzmann_density = boole_write_boltzmann_density
    in%boole_write_electric_potential = boole_write_electric_potential
    in%boole_write_moments = boole_write_moments
    in%boole_write_fourier_moments = boole_write_fourier_moments
    in%boole_write_exit_data = boole_write_exit_data
    in%boole_write_grid_data = boole_write_grid_data
    in%boole_preserve_energy_and_momentum_during_collisions = boole_preserve_energy_and_momentum_during_collisions
    in%n_electric_potential_updates = n_electric_potential_updates
    in%update_dimension = update_dimension
    in%n_species = n_species
    in%boole_static_ne = boole_static_ne

    print *,'GORILLA_APPLETS: Loaded input data from anomalous_transport.inp'

end subroutine read_anomalous_transport_inp_into_type

end module utils_anomalous_transport_mod
