%% 2d color plot
clear all
folder = '/temp/schatz_j/GORILLA_self_consistent_electric_field/GORILLA_APPLETS/EXAMPLES/self_consistent_electric_field/';

%Load metadata about grid
netcdf_file = [folder,'grid_data.nc'];
grid_size = ncread(netcdf_file,"grid_size");
ns = grid_size(1);
ntheta = grid_size(2);
nphi = grid_size(3);
n_extra_rings =  ncread(netcdf_file,"number of extra rings");

%Load grid data
coordinates = load([folder,'vertex_coordinates.dat']);
indices = load([folder,'vertex_indices.dat']);
volumes = load([folder,'prism_volumes.dat']);
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
    mesh_y(1:3,t) = coordinates(indices_plot(t,:),2);
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
prism_moments = load([folder,'1_prism_moments.dat']);
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
n_prisms_inner_ring = n_extra_rings*ntheta*2;
volume_times_density_inner_ring = sum(moment_plot(1:n_prisms_inner_ring).*...
                                      volumes_plot(1:n_prisms_inner_ring));
volume_inner_ring = sum(volumes_plot(1:n_prisms_inner_ring));
density_inner_ring = volume_times_density_inner_ring / volume_inner_ring;
moment_plot(1:n_prisms_inner_ring) = density_inner_ring;

figure
fill(mesh_x,mesh_y,moment_plot,'LineWidth',mesh_thickness)
xlabel('$s$ / cm','interpreter','latex')
ylabel('$\theta$ / cm','interpreter','latex')
title('density')
grid on
colorbar
clim([0,5e13])

