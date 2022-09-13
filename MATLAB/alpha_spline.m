%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% Calculate one orbit with Gorilla for a passing Deuterium particle with a field-aligned grid for non-axisymmetric VMEC
% Compare the results for different polynominal orders in the calculation of Gorilla
% Create a figure with the poincare plot and the total energy
%
%#######################################################################################################################
% programed functions:
%-----------------------------------------------------------------------------------------------------------------------
% InputFile
% NameList
% read_in
% load_copy
%#######################################################################################################################
%author: Michael Scheidt
%created: 19.01.2021


%Initialize used paths and files
%-----------------------------------------------------------------------------------------------------------------------

%Name of the current calculation to create folder
name_test_case='alpha_spline';

%path of Matlab script
path_script=pwd;

%main path of Gorilla
c=strsplit(path_script,'/');
path_main=strjoin(c(1:end-1),'/');

%path to Run Gorilla code in
mkdir(path_main,['RUN/',name_test_case]);
path_RUN=[path_main,'/RUN/',name_test_case];

%path of blueprints of input files
path_inp_blueprints=[path_main,'/inp_blueprints'];

%path of the used functions
path_functions=[path_main,'/MATLAB/functions'];

%define path for data and plots and create a new folder their
mkdir([path_script,'/data_plots'],name_test_case);
path_data_plots=[path_script,'/data_plots/',name_test_case];

%Change to RUN folder and clean it
cd(path_RUN);

%add path
addpath(path_inp_blueprints);
addpath(path_functions);

%initialize class for inputfiles
gorilla = InputFile([path_inp_blueprints,'/gorilla.inp']);
gorilla_plot = InputFile([path_inp_blueprints,'/gorilla_plot.inp']);
tetra_grid = InputFile([path_inp_blueprints,'/tetra_grid.inp']);

%read default inputfiles
%All variables get a default value from the blueprints
gorilla.read();
gorilla_plot.read();
tetra_grid.read();


%Set input variables for Gorilla
%-----------------------------------------------------------------------------------------------------------------------

%Input file gorilla
%All necessary variables for current calculation

    %Change in the electrostatic potential within the plasma volume in Gaussian units
        gorilla.GORILLANML.eps_Phi=0;

    %Coordinate system
        %1 ... (R,phi,Z) cylindrical coordinate system
        %2 ... (s,theta,phi) symmetry flux coordinate system
        gorilla.GORILLANML.coord_system = 2;

    %particle species
        %1 ... electron, 2 ... deuterium ion, 3 ... alpha particle
        gorilla.GORILLANML.ispecies = 3;

    %Switch for initial periodic coordinate particle re-location (modulo operation)
        % true ... Particles are re-located at initialization in the case of a periodic coordinate, if they are outside the computation domain.
        % false ... Particles are not re-located at initialization (This might lead to error if particles are outside the computation domain)
        gorilla.GORILLANML.boole_periodic_relocation = true;

    %1 ... numerical RK pusher, 2 ... polynomial pusher
        gorilla.GORILLANML.ipusher = 2;

    %Polynomial order for orbit pusher (from 2 to 4)
        gorilla.GORILLANML.poly_order = 4;


%Input file tetra_grid
%All necessary variables for current calculation

    %Grid Size
        %Rectangular: nR, Field-aligned: ns
        tetra_grid.TETRA_GRID_NML.n1 = 100;
        %Rectangular: nphi, Field-aligned: nphi
        tetra_grid.TETRA_GRID_NML.n2 = 100;
        %Rectangular: nZ, Field-aligned: ntheta
        tetra_grid.TETRA_GRID_NML.n3 = 100;

    %Grid kind
        %1 ... rectangular grid for axisymmetric EFIT data
        %2 ... field-aligned grid for axisymmetric EFIT data
        %3 ... field-aligned grid for non-axisymmetric VMEC
        tetra_grid.TETRA_GRID_NML.grid_kind = 3;

    %Switch for selecting number of field periods automatically or manually
        %.true. ... number of field periods is selected automatically (Tokamak = 1, Stellarator depending on VMEC equilibrium)
        %.false. ... number of field periods is selected manually (see below)
        tetra_grid.TETRA_GRID_NML.boole_n_field_periods = true;


%Input file gorilla_plot
%All necessary variables for current calculation

    %Switch for options
        % 1 ... Single orbit - Starting positions and pitch for the orbit are taken from file (see below) [First Line]
        % 2 ... Single orbit - Starting positions and pitch for the orbit are taken from starting drift surfaces (see below)
        % 3 ... Multiple orbits - Starting positions and pitch for orbits are taken from file (see below) [Every Line New Starting position]
        % 4 ... Multiple orbits - Starting positions and pitch for orbits are taken from drift surfaces with regular spacing (see below)
        gorilla_plot.GORILLA_PLOT_NML.i_orbit_options = 2;

    %Total individual orbit flight time for plotting
        gorilla_plot.GORILLA_PLOT_NML.total_orbit_time = 0.0008;

    %Total Energy of particle in eV
        gorilla_plot.GORILLA_PLOT_NML.energy_eV_start = 3.5.*10.^6;

    %Switch for plotting Poincar� cuts at toroidal variable $\varphi$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_phi_0 = true;

        %Number of skipped (non-printed) Poincar� cuts at toroidal variable $\varphi$ = 0
            gorilla_plot.GORILLA_PLOT_NML.n_skip_phi_0 = 1;

        %Filename for Poincar� cuts at toroidal variable $\varphi$ = 0 in cylindrical coordinates (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_rphiz = 'poincare_plot_phi_0_rphiz.dat';

        %Filename for Poincar� cuts at toroidal variable $\varphi$ = 0 in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
            gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_sthetaphi = 'poincare_plot_phi_0_sthetaphi.dat';

    %Switch for plotting Poincar� cuts at parallel velocity $v_\parallel$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_vpar_0 = false;

    %Switch for plotting full orbit
        gorilla_plot.GORILLA_PLOT_NML.boole_full_orbit = false;

        %Number of skipped (non-printed tetrahedra passings) full orbit
            gorilla_plot.GORILLA_PLOT_NML.n_skip_full_orbit = 1;

        %Filename for full orbit in cylindrical coordinates (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_rphiz = 'full_orbit_plot_rphiz.dat';

        %Filename for full orbit in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
            gorilla_plot.GORILLA_PLOT_NML.filename_full_orbit_sthetaphi = 'full_orbit_plot_sthetaphi.dat';

    %Plot invariances of motion (ONLY for single orbits)

        %Switch for plotting total particle energy
            gorilla_plot.GORILLA_PLOT_NML.boole_e_tot = false;

        %Switch for plotting canoncial (toroidal) angular momentum $p_\varphi$
            gorilla_plot.GORILLA_PLOT_NML.boole_p_phi = false;

        %Switch for parallel adiabatic invariant $J_\parallel$
            gorilla_plot.GORILLA_PLOT_NML.boole_J_par = false;

    %Single orbit from starting drift surface (i_orbit_options = 2)

        %Starting drift surface
            % = s for (s,$\vartheta$,$\varphi$)
            % = R for (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.start_pos_x1_beg = 0.5;

        %End drift surface
            % = s for (s,$\vartheta$,$\varphi$)
            % = R for (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.start_pos_x1_end = 0.9;

        %Number of drift surfaces in between start and end
            gorilla_plot.GORILLA_PLOT_NML.n_surfaces = 1;

        %Starting value for toroidal variable
            % = $\vartheta$ for (s,$\vartheta$,$\varphi$)
            % = $\varphi$ for (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.start_pos_x2  = 0.1;

        %Starting value for poloidal variable $\vartheta$
            % = $\varphi$ for (s,$\vartheta$,$\varphi$)
            % = Z for (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.start_pos_x3  = 0.1;

        %Pitch parameter $\lambda$ = $v_\parallel$ / vmod
            gorilla_plot.GORILLA_PLOT_NML.start_pitch_parameter = 0.95;


%Run Gorilla
%-----------------------------------------------------------------------------------------------------------------------

%write Input files for Gorilla
gorilla.write([path_RUN,'/gorilla.inp']);
gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);
tetra_grid.write([path_RUN,'/tetra_grid.inp']);

%Create softlinks for used files
! ln -s ../../test_gorilla_main.x .
! ln -s ../../DATA .
! ln -s DATA/VMEC/netcdf_file_for_test.nc .
! ln -s DATA/VMEC/netcdf_file_for_test_6.nc .

%Necessary files for Tokamak
! ln -s ../../inp_blueprints/field_divB0.inp .
! ln -s ../../inp_blueprints/preload_for_SYNCH.inp .
! ln -s ../../inp_blueprints/seed.inp .

%Run Gorilla Code
! ./test_gorilla_main.x

data=struct();
data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);
%Second Run for trapped particle
%-----------------------------------------------------------------------------------------------------------------------



%Create Plots of created data
%-----------------------------------------------------------------------------------------------------------------------

theta_spline=linspace(min(data.poincare_plot_phi_0_sthetaphi(:,2)),max(data.poincare_plot_phi_0_sthetaphi(:,2)),100.*numel(data.poincare_plot_phi_0_sthetaphi(:,2)));
s_spline=spline(data.poincare_plot_phi_0_sthetaphi(:,2),data.poincare_plot_phi_0_sthetaphi(:,1),theta_spline);

gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_rphiz = 'poincare_plot_phi_0_rphiz_long.dat';
gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_sthetaphi = 'poincare_plot_phi_0_sthetaphi_long.dat';
gorilla_plot.GORILLA_PLOT_NML.total_orbit_time = 1;

gorilla.write([path_RUN,'/gorilla.inp']);
gorilla_plot.write([path_RUN,'/gorilla_plot.inp']);
tetra_grid.write([path_RUN,'/tetra_grid.inp']);
! ./test_gorilla_main.x

data=struct();
data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);

figure
plot(data.poincare_plot_phi_0_sthetaphi_long(:,2),data.poincare_plot_phi_0_sthetaphi_long(:,1),'s','MarkerSize',4)
hold on
plot(theta_spline,s_spline)
xlabel('\vartheta')
ylabel('s')
title('Poincar\.{e} plot ($\varphi=0$) Stellerator','Interpreter','latex');
hold off

savefig([path_data_plots,'/poincare_plot.fig']);
saveas(gcf,[path_data_plots,'/poincare_plot.png']);


%Go back to path of Matlab script
cd(path_script);




