%% 2d color plot
clear all
folder = '/temp/schatz_j/GORILLA_self_consistent_electric_field/GORILLA_APPLETS/EXAMPLES/self_consistent_electric_field/';
folder = [pwd,'/'];

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

%% 1d plot
clc

folder = '/proj/plasma/CODE/schatz_j/GORILLA_APPLETS/EXAMPLES/self_consistent_electric_field/';
%folder = '/temp/schatz_j/self_consistent_electric_field_runs/flat_source/2025.09.13_1000_ions_small_factor_from_tracing_time/';
%folder = '/temp/schatz_j/self_consistent_electric_field_runs/flat_source/2025.09.15_1000_ions_small_factor_from_tracing_time_positive_potenital_at_boundary/';
folder = '/temp/schatz_j/self_consistent_electric_field_runs/flat_source/2025.09.15_1000_ions_small_factor_from_tracing_time_positive_potenital_at_boundary_optimised_runtime/';
%folder = '/temp/schatz_j/self_consistent_electric_field_runs/flat_source/2025.09.16_1000_ions_small_factor_from_tracing_time_ultra_positive_potenital_at_boundary_optimised_runtime/';
%folder = '/temp/schatz_j/self_consistent_electric_field_runs/flat_source/2025.09.13_1000_ions_large_factor_from_tracing_time/';
%folder = '/itp/MooseFS/schatz_j/self_consistent_electric_field_runs/center_source/2025.09.08/';
%folder = '/temp/schatz_j/self_consistent_electric_field_runs/honest_electrons_central_source (electron density scaled wrongly, background density linear in s))/';
%folder = '/temp/schatz_j/self_consistent_electric_field_runs/flat_source/2025.09.14_1000_ions_small_factor_from_tracing_time_zero_potenital_at_boundary/';
%folder = '/temp/schatz_j/self_consistent_electric_field_runs/frozen_electrons/2025.09.09_1000_ions/';
%folder = '/temp/schatz_j/self_consistent_electric_field_runs/frozen_electrons/2025.09.16_1000_ions_used_for_poster/';
folder = [pwd,'/'];


volumes = load([folder,'prism_volumes.dat']);

N = 10;
N_low = 2;

ns = 30;
nphi = 30;
ntheta = 30;
n_triangles = 30*30*2;

one_d_plot = zeros(ns,N-N_low+1);
potential_plot =  zeros(ns+1,N-N_low+1);
f_b_plot = zeros(ns+1,N-N_low+1);
echarge = 4.8032e-10;
T = 3.5e3*1.6022e-12;

s_shell_volumes = zeros(ns,1);
for i = 1:ns
    for j = 1:nphi
        range_low  = (j-1)*ntheta*ns*2+(i-1)*2*ntheta+1;
        range_high = (j-1)*ntheta*ns*2+ i   *2*ntheta;
        s_shell_volumes(i) = s_shell_volumes(i) + sum(volumes(range_low:range_high));
    end
end

for j = 1:N-N_low+1
    density = load([folder,'ion_densities_after_electric_potential_update_',num2str(N_low+j-1),'.dat']);
    num_ions_plot = zeros(n_triangles,1);
    for i = 1:nphi
        range = 1+(i-1)*n_triangles:i*n_triangles;
        num_ions_plot = num_ions_plot + density(range).*volumes(range);
    end
    for k = 1:ns
        range = (k-1)*ntheta*2+1:k*ntheta*2;
        one_d_plot(k,j) = sum(num_ions_plot(range))/s_shell_volumes(k);
        %one_d_plot(k,j) = sum(num_ions_plot((k-1)*ntheta*2+1:k*ntheta*2))/(2*ntheta);
    end
    potential = load([folder,'phi_elec_after_electric_potential_update_',num2str(N_low+j-1),'.dat']);
    potential = potential(1:30:30*31);%+potential(31:30:30*31))/2;
    f_b_plot(:,j) = exp(-echarge*potential/T)*3e13;
    potential_plot(:,j) = potential;
end

starting_potential =  load([folder,'phi_elec_after_electric_potential_update_1.dat']);
start_potential_plot = starting_potential(1:30:31*30);%+starting_potential(31:30:30*31))/2;
starting_distribution = load([folder,'ion_densities_after_electric_potential_update_1.dat']);
start_2d = zeros(n_triangles,1);
start_plot = zeros(ns,1);
for i = 1:nphi
    range = 1+(i-1)*n_triangles:i*n_triangles;
    start_2d = start_2d + starting_distribution(range).*volumes(range);
end
for k = 1:ns
    range = (k-1)*ntheta*2+1:k*ntheta*2;
    start_plot(k) = sum(start_2d(range))/s_shell_volumes(k);
    %start_plot(k) = sum(start_2d((k-1)*ntheta*2+1:k*ntheta*2))/(2*ntheta);
end

s = linspace(1e-4,1,31);
s_center = (s(2:end)+s(1:end-1))/2;
electron_density_factor = 1.8129425031187620;
electron_density = (1-0.9*s_center)*1e13*electron_density_factor;
electron_density = electron_density';

%figure
%plot(start_plot(:),'r')
%hold on
%plot(1,0)
%for i = 1:N-N_low+1
%    plot(one_d_plot(:,i),'b')
%    %plot(f_b_plot(:,i),'b:')
%end
%%%plot(electron_density,'m')
%hold off
%grid on
%ylabel('ion density')

%average_electron_density = sum(electron_density.*s_shell_volumes)/sum(s_shell_volumes)
%for i = 1:N-N_low+1
%average_ion_density = sum(one_d_plot(:,i).*s_shell_volumes)/sum(s_shell_volumes)
%end

figure
plot(start_potential_plot,'r')
hold on
plot(1,0)
for i = 1:N-N_low+1
    plot(potential_plot(:,i),'b')
end
hold off
grid on
ylabel('potential')
%ylim([-60,40])


one_d = load([folder,'one_d_densities1.dat']);
one_d_save = one_d;
figure
plot(one_d(:,2),'b')
hold on
plot(one_d(:,1),'r')
 for i = N_low:N
     one_d = load([folder,'one_d_densities',num2str(i),'.dat']);
     plot(one_d(:,1),'m')
     plot(one_d(:,2),'k')
 end
%plot(f_b_plot(:,N-N_low+1),'r')
hold off
grid on

 %% quadratic fit for diffusion coefficient
folder = [pwd,'/'];
data = load([folder,'A_and_B.dat']); 

B = data(:,2);
x = linspace(1e-4,1,31);
x_extend = linspace(0,1);
nstart = 1;
nend = 26;
B=B(nstart:nend);
x=x(nstart:nend);

p = polyfit(x, B, 2)

coeff = [2.2745400112327060E-002,  0.12463948389456188,       0.15247066690352659, 0];

a = figure;
plot(x,B,'x')
hold on
plot(x_extend,coeff(1)+coeff(2)*x_extend+coeff(3)*x_extend.^2)
%plot(x,p(3)+p(2)*x+p(1)*x.^2)
hold off
grid on
xlabel('flux surface label $s$', 'interpreter','latex')
ylabel('$D_{11}^{ss}$', 'interpreter','latex')
set(a,'Units','centimeters','position',[20,10,12,8.8]);

name = 'diffusion_coefficient_fit.png';


%% plot s-values over time

folder = [pwd,'/'];
folder = '/temp/schatz_j/self_consistent_electric_field_runs/s-development_for_plots_in_dissertation/18.11.2025_long_tracing_time/';

istart = 15;
iend = 15;
n = 3;
for i = istart:iend
    data = load([folder,'data',num2str(i),'.dat']);
    tau = data(end-1,1);
    data(end,:) = [];
    data(:,1)= data(:,1)/tau;
    xfit = data(numel(data(:,1))*0.5:end,1);
    yfit = data(numel(data(:,1))*0.5:end,n);
    p = polyfit(xfit, yfit, 1);
    a = figure
    plot(data(:,1),data(:,n))
    hold on
    plot(data(:,1),p(2)+p(1)*data(:,1),'r')
    hold off
    grid on
    ylabel('$\langle(s(t)-s_0)^2\rangle$','interpreter','Latex')
    xlabel('$t / \tau_{\mathrm{tr}}$','interpreter','Latex')
    if n==2
        ylabel('$\langle s(t)-s_0 \rangle$','interpreter','Latex')
    elseif n==3
        ylabel('$\langle(s(t)-s_0)^2\rangle$','interpreter','Latex')
    end
    %xlim([0,data(end-1,1)])
    set(a,'Units','centimeters','position',[20,10,8,8.8]);
end

if n==2
    name = 'deltas_mid_radius.png';
elseif n==3
    name = 'deltas2_mid_radius.png';
end

%% plot s-values at t_end

folder = [pwd,'/'];
data = load('/temp/schatz_j/self_consistent_electric_field_runs/s-development_for_plots_in_dissertation/18.11.2025_long_tracing_time/exit_s_values15.dat');

nbins = 50;
nx_analytical = 10000;
nratio = nx_analytical/nbins;
srange = max(data)-min(data);
smin  = min(data)-srange*0.05;
smax  = max(data)+srange*0.05;
edges = linspace(smin, smax, nbins + 1);  % Nbins intervals → Nbins+1 edges
centers = (edges(1:end-1)+edges(2:end))/2;
s0 = 0.5; % your desired center
[~, idx] = min(abs(centers - s0)); % Find the bin whose center is closest to s0
shift = s0 - centers(idx); % Compute the shift needed to move that bin center to s0
edges = edges + shift;
%histogram(x, 'BinEdges', edges);

a = figure
h = histogram(data,'BinEdges', edges);%,'Edgecolor','none');
set(gca, 'YScale', 'log');
hold on
%set(gca,'YScale','log') 
mean = sum(data)/numel(data)
for i = 1:numel(data)
    if abs(data(i)-mean)>0.002
        %data(i)=mean;
    end
end

ordered_data = sort(abs(data-mean));
standard_deviation1 = ordered_data(floor(numel(data)*0.6827))
standard_deviation2 = sqrt(sum(abs(data-mean).^2)/numel(data))
x = linspace(min(data),max(data),nx_analytical);
y = exp(-(x-mean).^2/(2*standard_deviation2^2));
normalisation = sum(y)/numel(data)/nratio;
y = y/normalisation;
semilogy(x,y,'linewidth',1)
hold off
ylim([1e-1,1e5])
xlabel('flux surface label s')
ylabel('number of electrons / 1')
grid on
ax = gca;
ax.XMinorGrid = 'off';    % turn off minor grid on x-axis
ax.YMinorGrid = 'off';    % turn off minor grid on y-axis
set(a,'Units','centimeters','position',[20,10,18,8.8]);

j = 0;
for i = 1:numel(data)
    if abs(data(i)-mean)>0.005
       j = j+1;
    end
end
result = j;

name = 'final_s_distribution_mid_radius.png';

%% single particle s-values over time

single_particle_data = load('/afs/itp.tugraz.at/proj/plasma/CODE/schatz_j/GORILLA_APPLETS/EXAMPLES/self_consistent_electric_field/fort.71');

figure
plot(single_particle_data(:,1),single_particle_data(:,2))
xlabel('t')
ylabel('s')
grid on

figure
plot(single_particle_data(:,1),single_particle_data(:,3),'r')
hold on
plot(single_particle_data(:,1),single_particle_data(:,4),'b')
xlabel('t')
ylabel('theta, phi')
grid on

figure
plot(single_particle_data(:,1),single_particle_data(:,5),'r')
hold on
plot(single_particle_data(:,1),single_particle_data(:,6),'b')
%plot(single_particle_data(:,1),single_particle_data(:,6).^2+single_particle_data(:,5).^2,'m')
hold off
xlabel('t')
ylabel('vpar, vperp,vmod')
grid on


%% Is this the data from ALEX RUNOV?
rad = linspace(0,1);
%Electron density [/10^20/m^3]:
  dens_E20 = 0.37d0*(0.74d0 + 0.26d0*(1.d0-rad.^2.5).^1.5 - 0.06d0*(1.d0-exp(-rad.^2*4.d0)));
  density = dens_E20*1.d14;  % cm^{-3}
%Ion temperature [keV]:
  ti_keV = 1.77d0*(0.06d0 + 0.94d0*(1.d0 - rad.^2.9).^2);
  ti_eV = ti_keV*1.d3;
% Electron temperature [keV]:
  te_keV = ti_keV + 3.19d0*(1.d0 - rad.^2.1).^5.4;
  te_eV = te_keV*1.d3;

  figure
  plot(rad, density)
  title('density')
  grid on

  figure
  plot(rad,ti_keV)
  title('ion temperature keV')
  grid on

  figure
  plot(rad,te_keV)
  title('electron temperature keV')
  grid on