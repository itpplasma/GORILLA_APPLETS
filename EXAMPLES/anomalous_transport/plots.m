%% 2d color plot
clear all

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

%'x_prism_moments.dat' contains prism moments after the last potential
%update, for ions x=1, for electrons x=2
%Currently, electrons are computed with random walk, therefore only
%electron densities are recorded, other moments are set to zero
prism_moments = load('prism_moments.dat');
%All odd colums conatain real data, all even columns the corresponding
%imaginary data, for the current settings all even columns are zero
%the switches in "gorilla.inp" to turn on moment computation are called
%boole_time_hamiltonian, boole_gyrophase, boole_vpar_int, boole_vpar2_int
%if all are set to .true., n=1 selects density and n=5 selects parallel 
%velocity (not normalised by density and charge)
n=1;
moment = prism_moments(:,n);
moment_plot = zeros(n_triangles,1);
for i = 1:nphi
    moment_plot = moment_plot + moment(1+(i-1)*n_triangles:i*n_triangles)/nphi;
end
%for inner rings with very small volumes, perform an average
n_prisms_inner_ring = n_extra_rings*nz*2;
volume_times_density_inner_ring = sum(moment_plot(1:n_prisms_inner_ring).*...
                                      volumes_plot(1:n_prisms_inner_ring));
volume_inner_ring = sum(volumes_plot(1:n_prisms_inner_ring));
density_inner_ring = volume_times_density_inner_ring / volume_inner_ring;
moment_plot(1:n_prisms_inner_ring) = density_inner_ring;

exit_data = load('exit_data.dat');

figure
fill(mesh_x,mesh_y,moment_plot,'LineWidth',mesh_thickness)
% hold on
% plot(exit_data(:,4),exit_data(:,6), 'rx')
% hold off
xlabel('$R$ / cm','interpreter','latex')
ylabel('$Z$ / cm','interpreter','latex')
title('density')
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
title('1D Radial Density Profile (flux surface averaged)')
grid on
ylim([0,5]*1e13)

%% Diffusion coefficient analysis
% Load diffusion coefficient data
diffusion_data = load('diffusion_coefficient_data.dat');

% Extract columns: time, <delta_s>, <delta_s^2>, n_samples
time = diffusion_data(:,1);
delta_s = diffusion_data(:,2);
delta_s_squared = diffusion_data(:,3);
n_samples = diffusion_data(:,4);

% Plot displacement vs time
figure
subplot(2,1,1)
plot(time*1e6, delta_s, 'b-', 'LineWidth', 2)
xlabel('Time ($\mu$s)', 'interpreter', 'latex')
ylabel('$\langle s - s_0 \rangle$', 'interpreter', 'latex')
title('Mean Radial Displacement vs Time')
grid on

% Plot squared displacement vs time
subplot(2,1,2)
plot(time*1e6, delta_s_squared, 'r-', 'LineWidth', 2)
hold on
% Add linear fit to show diffusive behavior: <delta_s^2> = 2*B*t
p = polyfit(time, delta_s_squared, 1);
plot(time*1e6, polyval(p, time), 'k--', 'LineWidth', 1.5)
xlabel('Time ($\mu$s)', 'interpreter', 'latex')
ylabel('$\langle (s - s_0)^2 \rangle$', 'interpreter', 'latex')
title('Mean Squared Radial Displacement vs Time')
legend('Data', sprintf('Linear fit: B = %.2e s^2/s', p(1)/2), 'Location', 'southeast')
grid on

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
title('Distribution of Particles on Flux Surfaces at Exit')
grid on
xlim([0 1])

% Add statistics to the plot
mean_s = mean(valid_flux);
std_s = std(valid_flux);
text(0.6, 0.9, sprintf('Mean: %.3f\nStd: %.3f', mean_s, std_s), ...
     'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w')

