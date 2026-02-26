%% Helical Core Particle Tracing - Visualization
% This script visualizes results from helical core particle tracing.
% It is a simplified version of anomalous_transport/plots.m without
% diffusion coefficient analysis or electric potential perturbation plots.

clear all

%% User settings
% Select which quantity to plot in the 2D color plot:
%   'density'  - Particle density (t_hamiltonian, moment 1)
%   'vpar'     - Parallel velocity (vpar_int, moment 3)
%   'vpar2'    - Parallel velocity squared (vpar2_int, moment 4)
plot_quantity = 'vpar';

%% 2D color plot of selected quantity
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

%Load prism moments
% Column structure (with all 4 optional quantities enabled in gorilla.inp):
%   Columns 1-2: t_hamiltonian (real, imag) - proportional to density
%   Columns 3-4: gyrophase (real, imag)
%   Columns 5-6: vpar_int (real, imag) - parallel velocity integral
%   Columns 7-8: vpar2_int (real, imag) - parallel velocity squared integral
prism_moments = load('prism_moments.dat');

% Select moment column based on plot_quantity setting
switch plot_quantity
    case 'density'
        moment_col = 1;  % t_hamiltonian (real part)
        plot_title = 'Particle Density - Helical Core';
        colorbar_label = 'Density';
        clim_values = [0, 5e13];
    case 'vpar'
        moment_col = 5;  % vpar_int (real part)
        plot_title = 'Parallel Velocity Distribution - Helical Core';
        colorbar_label = '$\langle v_\parallel \rangle$';
        clim_values = 'auto';
    case 'vpar2'
        moment_col = 7;  % vpar2_int (real part)
        plot_title = 'Parallel Velocity Squared - Helical Core';
        colorbar_label = '$\langle v_\parallel^2 \rangle$';
        clim_values = 'auto';
    otherwise
        error('Unknown plot_quantity: %s. Use ''density'', ''vpar'', or ''vpar2''.', plot_quantity);
end

moment = prism_moments(:, moment_col);
moment_plot = zeros(n_triangles,1);
for i = 1:nphi
    moment_plot = moment_plot + moment(1+(i-1)*n_triangles:i*n_triangles)/nphi;
end
%For inner rings with very small volumes, perform an average
n_prisms_inner_ring = n_extra_rings*nz*2;
volume_times_moment_inner_ring = sum(moment_plot(1:n_prisms_inner_ring).*...
                                      volumes_plot(1:n_prisms_inner_ring));
volume_inner_ring = sum(volumes_plot(1:n_prisms_inner_ring));
moment_inner_ring = volume_times_moment_inner_ring / volume_inner_ring;
moment_plot(1:n_prisms_inner_ring) = moment_inner_ring;

exit_data = load('exit_data.dat');

figure
fill(mesh_x,mesh_y,moment_plot)%,'EdgeColor','none')
xlabel('$R$ / cm','interpreter','latex')
ylabel('$Z$ / cm','interpreter','latex')
title(plot_title)
grid on
cb = colorbar;
cb.Label.String = colorbar_label;
cb.Label.Interpreter = 'latex';
if ~strcmp(clim_values, 'auto')
    clim(clim_values);
end

%% 1D radial profile (volume-averaged over flux surfaces)
% Number of triangles per radial shell (2 triangles per quadrilateral cell in Z)
n_triangles_per_radial = nz * 2;

% Compute volume-weighted average for each radial shell
moment_radial = zeros(nr, 1);
volume_radial = zeros(nr, 1);

for ir = 1:nr
    % Triangle indices for this radial shell
    t_start = (ir - 1) * n_triangles_per_radial + 1;
    t_end = ir * n_triangles_per_radial;

    % Accumulate volume-weighted moment and total volume
    moment_radial(ir) = sum(moment_plot(t_start:t_end) .* volumes_plot(t_start:t_end));
    volume_radial(ir) = sum(volumes_plot(t_start:t_end));
end

% Normalize by volume to get average in each shell
moment_radial_avg = moment_radial ./ volume_radial;

% Set axis labels based on plot_quantity
switch plot_quantity
    case 'density'
        radial_ylabel = 'Volume-averaged density';
        radial_title = '1D Radial Density Profile - Helical Core';
        radial_ylim = [0, 5]*1e13;
    case 'vpar'
        radial_ylabel = 'Volume-averaged $\langle v_\parallel \rangle$';
        radial_title = '1D Radial Parallel Velocity Profile - Helical Core';
        radial_ylim = 'auto';
    case 'vpar2'
        radial_ylabel = 'Volume-averaged $\langle v_\parallel^2 \rangle$';
        radial_title = '1D Radial $v_\parallel^2$ Profile - Helical Core';
        radial_ylim = 'auto';
end

figure
plot(1:nr, moment_radial_avg, 'b-', 'LineWidth', 2)
xlabel('Radial index', 'interpreter', 'latex')
ylabel(radial_ylabel, 'interpreter', 'latex')
title(radial_title)
grid on
if ~strcmp(radial_ylim, 'auto')
    ylim(radial_ylim);
end

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
