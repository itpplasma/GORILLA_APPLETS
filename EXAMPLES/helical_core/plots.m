%% Helical Core Particle Tracing - Visualization
% This script visualizes results from helical core particle tracing.
% It is a simplified version of anomalous_transport/plots.m without
% diffusion coefficient analysis or electric potential perturbation plots.

clear all

%% 2D color plot of density
%Load metadata about grid
netcdf_file = 'grid_data.nc';
grid_size = ncread(netcdf_file,"grid_size");
nr = grid_size(1);
nphi = grid_size(2);
nz = grid_size(3);
n_extra_rings =  ncread(netcdf_file,"number of extra rings");

%Load grid data
coordinates = load('vertex_coordinates.dat');
indices = load('vertex_indices.dat');
volumes = load('prism_volumes.dat');
nvert = numel(coordinates(:,1));
ntetr = numel(indices(:,1));

%Compute mesh for plotting
phi = 0.5; %determines at what angle of phi a cross-section is taken
ntetr_per_slice = ntetr/nphi;
n_triangles = ntetr_per_slice/3;
slice = floor(phi/((2*pi/5)/nphi))+1;
tetr_plot = ntetr_per_slice*(slice-1)+1:3:ntetr_per_slice*slice;
indices_plot = indices(tetr_plot,1:3);
volumes_plot = volumes(n_triangles*(slice-1)+1:n_triangles*slice);
mesh_x = zeros(4,n_triangles);
mesh_y = zeros(4,n_triangles);
for t = 1:n_triangles
    mesh_x(1:3,t) = coordinates(indices_plot(t,:),1);
    mesh_y(1:3,t) = coordinates(indices_plot(t,:),3);
    mesh_x(4,t) = mesh_x(1,t);
    mesh_y(4,t) = mesh_y(1,t);

    if (max(mesh_y(:,t)) - min(mesh_y(:,t))) > pi
        for i = 1:4
            if mesh_y(i,t)==0
                mesh_y(i,t)=2*pi;
            end
        end
    end
end
mesh_color = [204,204,204]/256;
mesh_thickness = 0.01;

%Load prism moments (density distribution)
prism_moments = load('prism_moments.dat');
%Column 1 contains real density data (column 2 is imaginary, should be zero)
n=1;
moment = prism_moments(:,n);
moment_plot = zeros(n_triangles,1);
for i = 1:nphi
    moment_plot = moment_plot + moment(1+(i-1)*n_triangles:i*n_triangles)/nphi;
end
%For inner rings with very small volumes, perform an average
n_prisms_inner_ring = n_extra_rings*nz*2;
volume_times_density_inner_ring = sum(moment_plot(1:n_prisms_inner_ring).*...
                                      volumes_plot(1:n_prisms_inner_ring));
volume_inner_ring = sum(volumes_plot(1:n_prisms_inner_ring));
density_inner_ring = volume_times_density_inner_ring / volume_inner_ring;
moment_plot(1:n_prisms_inner_ring) = density_inner_ring;

exit_data = load('exit_data.dat');

figure
fill(mesh_x,mesh_y,moment_plot,'LineWidth',mesh_thickness)
xlabel('$R$ / cm','interpreter','latex')
ylabel('$Z$ / cm','interpreter','latex')
title('Particle Density - Helical Core')
grid on
colorbar
clim([0,5e13])

%% 1D radial density profile (volume-averaged over flux surfaces)
% Number of triangles per radial shell (2 triangles per quadrilateral cell in Z)
n_triangles_per_radial = nz * 2;

% Compute volume-weighted average density for each radial shell
density_radial = zeros(nr, 1);
volume_radial = zeros(nr, 1);

for ir = 1:nr
    % Triangle indices for this radial shell
    t_start = (ir - 1) * n_triangles_per_radial + 1;
    t_end = ir * n_triangles_per_radial;

    % Accumulate volume-weighted density and total volume
    density_radial(ir) = sum(moment_plot(t_start:t_end) .* volumes_plot(t_start:t_end));
    volume_radial(ir) = sum(volumes_plot(t_start:t_end));
end

% Normalize by volume to get average density in each shell
density_radial_avg = density_radial ./ volume_radial;

figure
plot(1:nr, density_radial_avg, 'b-', 'LineWidth', 2)
xlabel('Radial index', 'interpreter', 'latex')
ylabel('Volume-averaged density', 'interpreter', 'latex')
title('1D Radial Density Profile - Helical Core')
grid on
ylim([0,5]*1e13)

%% Flux surface distribution at exit
% Load exit data
exit_data = load('exit_data.dat');

% Extract flux surface labels (last column)
flux_surface_labels = exit_data(:, end);

% Remove invalid entries (marked as -1 for particles outside domain)
valid_flux = flux_surface_labels(flux_surface_labels >= 0);

% Plot histogram
figure
histogram(valid_flux, 30, 'Normalization', 'probability', 'FaceColor', [0.2 0.6 0.8])
xlabel('Flux Surface Label $s$', 'interpreter', 'latex')
ylabel('Probability', 'interpreter', 'latex')
title('Distribution of Particles on Flux Surfaces at Exit - Helical Core')
grid on
xlim([0 1])

% Add statistics to the plot
mean_s = mean(valid_flux);
std_s = std(valid_flux);
text(0.6, 0.9, sprintf('Mean: %.3f\nStd: %.3f', mean_s, std_s), ...
     'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w')

%% Particle trajectories in R-Z plane
% Plot exit positions in the R-Z plane
figure
plot(exit_data(:,4), exit_data(:,6), 'b.', 'MarkerSize', 3)
xlabel('$R$ / cm', 'interpreter', 'latex')
ylabel('$Z$ / cm', 'interpreter', 'latex')
title('Final Particle Positions - Helical Core')
grid on
axis equal
