%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% Calculate one orbit with Gorilla for a passing Deuterium particle with a field-aligned grid for axisymmetric EFIT
% Compare the results for different polynominal orders in the calculation of Gorilla
% Create a figure with the Poincare plot and the total energy
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
%created: 19.02.2021


%Initialize used paths and files
%-----------------------------------------------------------------------------------------------------------------------

%Name of the current calculation to create folder
name_test_case='example_4';

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
        gorilla.GORILLANML.ispecies = 2;

    %Switch for initial periodic coordinate particle re-location (modulo operation)
        % true ... Particles are re-located at initialization in the case of a periodic coordinate, if they are outside the computation domain.
        % false ... Particles are not re-located at initialization (This might lead to error if particles are outside the computation domain)
        gorilla.GORILLANML.boole_periodic_relocation = false;

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
        tetra_grid.TETRA_GRID_NML.n2 = 30;
        %Rectangular: nZ, Field-aligned: ntheta
        tetra_grid.TETRA_GRID_NML.n3 = 30;

    %Grid kind
        %1 ... rectangular grid for axisymmetric EFIT data
        %2 ... field-aligned grid for axisymmetric EFIT data
        %3 ... field-aligned grid for non-axisymmetric VMEC
        tetra_grid.TETRA_GRID_NML.grid_kind = 2;

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
        gorilla_plot.GORILLA_PLOT_NML.i_orbit_options = 3;

    %Total individual orbit flight time for plotting
        gorilla_plot.GORILLA_PLOT_NML.total_orbit_time = 2;

    %Total Energy of particle in eV
        gorilla_plot.GORILLA_PLOT_NML.energy_eV_start = 3000;

    %Switch for plotting Poincar� cuts at toroidal variable $\varphi$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_phi_0 = true;

        %Number of skipped (non-printed) Poincar� cuts at toroidal variable $\varphi$ = 0
            gorilla_plot.GORILLA_PLOT_NML.n_skip_phi_0 = 100;

        %Filename for Poincar� cuts at toroidal variable $\varphi$ = 0 in cylindrical coordinates (R,$\varphi$,Z)
            gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_rphiz = 'poincare_plot_phi_0_rphiz.dat';

        %Filename for Poincar� cuts at toroidal variable $\varphi$ = 0 in symmetry flux coordinates (s,$\vartheta$,$\varphi$)
            gorilla_plot.GORILLA_PLOT_NML.filename_poincare_phi_0_sthetaphi = 'poincare_plot_phi_0_sthetaphi.dat';

    %Switch for plotting Poincar� cuts at parallel velocity $v_\parallel$ = 0
        gorilla_plot.GORILLA_PLOT_NML.boole_poincare_vpar_0 = false;

    %Switch for plotting full orbit
        gorilla_plot.GORILLA_PLOT_NML.boole_full_orbit = false;

    %Plot invariances of motion (ONLY for single orbits)

        %Switch for plotting total particle energy
            gorilla_plot.GORILLA_PLOT_NML.boole_e_tot = false;

        %Switch for plotting canoncial (toroidal) angular momentum $p_\varphi$
            gorilla_plot.GORILLA_PLOT_NML.boole_p_phi = false;

        %Switch for parallel adiabatic invariant $J_\parallel$
            gorilla_plot.GORILLA_PLOT_NML.boole_J_par = false;

    %Starting positions of particles (BOTH single and multiple orbits) and starting pitch parameter
        %File (either rphiz or sthetaphi) is automatically chosen dependent on coord_system in 'gorilla.inp')

        %Filename for list of starting position(s) of particle(s) in cylindrical coordinates (R,$\varphi$,Z) and pitch ($\lambda$)
            gorilla_plot.GORILLA_PLOT_NML.filename_orbit_start_pos_rphiz = 'orbit_start_rphizlambda.dat';

        %Filename for list of starting position(s) of particle(s) in symmetry flux coordinates (s,$\vartheta$,$\varphi$) and pitch ($\lambda$)
            gorilla_plot.GORILLA_PLOT_NML.filename_orbit_start_pos_sthetaphi = 'orbit_start_sthetaphilambda.dat';


%Create inputfile for starting positions
%s,theta,phi,lamda
start_pos_1=[0.8,0.0,0.63,0.7];
start_pos_2=[0.6,0.0,0.63,0.7];
start_pos_3=[0.4,0.0,0.63,0.4];
start_pos=[start_pos_1,start_pos_2,start_pos_3];
fileID = fopen([path_RUN,'/',gorilla_plot.GORILLA_PLOT_NML.filename_orbit_start_pos_sthetaphi],'w');
fprintf(fileID,'%.8f %.8f %.8f %.8f\n',start_pos);
fclose(fileID);

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

%Load Data-files and copy them into data_plots folder
data=struct();
data=load_copy(data,path_RUN,path_data_plots,gorilla_plot);

%Create Plots of created data
%-----------------------------------------------------------------------------------------------------------------------
%Poincare Plot passing particles in VMEC
figure
plot(data.poincare_plot_phi_0_rphiz(:,1),data.poincare_plot_phi_0_rphiz(:,3),'x')
hold on
xlabel('R [cm]')
ylabel('z [cm]')
title('Poincar\.{e} plot ($\varphi=0$) Tokamak','Interpreter','latex');
hold off

%Save figure as Matlab-file and as png
savefig([path_data_plots,'/poincare_plot.fig']);
saveas(gcf,[path_data_plots,'/poincare_plot.png']);

%Go back to path of Matlab script
cd(path_script);

