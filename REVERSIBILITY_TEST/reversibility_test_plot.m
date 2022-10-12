%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% * Compute collisionless guiding-center orbit with GORILLA for a trapped Deuterium particle with/without adaptive scheme.
% * Use a field-aligned grid for a non-axisymmetric VMEC MHD equilibrium.
% * Create a figure with the Poincar� plots (\varphi = 0) in cylindrical and symmetry flux coordinates.
% * Compute the normalized parallel adiabatic invariant as a function of banana bounces.
%
%#######################################################################################################################
% used functions:
%-----------------------------------------------------------------------------------------------------------------------
% InputFile
% NameList
% read_in
% load_copy
%#######################################################################################################################
% authors: Georg Graßler
% created: 16.09.2022

%For Construction of Demonstration
%-----------------------------------------------------------------------------------------------------------------------
%control of script
close all
boole_recalculate = true;
poly_order = 4;
grid_kind = 3;
[n1,n2,n3] = deal(70,70,70);
sfc_s_min = 1.d-3;
eps_Phi = -1.d-7;
%eps_Phi = 0;
ispecies = 2;
boole_guess = true;
i_time_tracing_option = 1;
t_total = 3.25d-5*6/7;
t_total = 1.3d-4*1/3.5;
dminor_0 = 0.25d0;
dminor_0 = 0.1d0;
phi_0 = 0.01;
contour_fraction = 1;
n_steps = 20;
n_orbits = 100;
n_snapshots = 2;
pitchpar_0 = 0.46d0;
pitchpar_0 = 0.6d0;
dpitchpar = 0.2d0;
energy_eV_0 = 3d3;
relative_bandwith = 0.1d0;
relative_bandwith = +0.1d0;
boole_diag_reversibility_test = false;

boole_show_contour = true;
boole_show_contour_back = false;
number_of_quantities = 7;
ylabel_quantities = {'$s$','$\vartheta$','$\varphi$','$\lambda$','$t$','$\phi$','$E$'};

n_skip_contour = 1;

forwardColor = [215,25,28]/256;
backwardColor = [44,123,182]/256;

MarkerContour = '-';
MarkerSizeContour = 15;
CentralLineWidth = 2;
FontSize = 25;
axis_factor = 3/5;

postprocessing_functions_path = 'postprocessing_functions';
addpath(postprocessing_functions_path);

%filenames for output
filename_reversibility_test = 'reversibility_test.dat';
filename_reversibility_test_back = 'reversibility_test_back.dat';

if (boole_recalculate)
  
    %Initialize used paths and files
    %-----------------------------------------------------------------------------------------------------------------------

    %Name of the current calculation to create folder
    name_test_case='reversibilty_test';

    %path of Matlab script
    path_script=pwd;

    %main path of GORILLA
    c=strsplit(path_script,'/');
    path_main=strjoin(c(1:end-1),'/');
    path_main_core = [path_main,'/../GORILLA'];

    %path to run GORILLA code
    mkdir(path_main,['EXAMPLES/MATLAB_RUN/',name_test_case]);
    path_RUN=[path_main,'/EXAMPLES/MATLAB_RUN/',name_test_case];

    %path of input files (blueprints)
    path_inp_files=[path_main,'/INPUT'];
    
    %path of input files of CORE GORILLA (blueprints)
    path_inp_files_core=[path_main_core,'/INPUT'];

    %path of the used functions
    path_functions=[path_main,'/MATLAB/functions'];

    %define path for data and plots and create a new folder
    mkdir([path_script,'/data_plots'],name_test_case);
    path_data_plots=[path_script,'/data_plots/',name_test_case];

    %Change to RUN folder and clean it
    cd(path_RUN);

    %add path
    addpath(path_inp_files);
    addpath(path_inp_files_core);
    addpath(path_functions);

    %initialize class for inputfiles
    gorilla = InputFile([path_inp_files_core,'/gorilla.inp']);
    tetra_grid = InputFile([path_inp_files_core,'/tetra_grid.inp']);
    gorilla_applets = InputFile([path_inp_files,'/gorilla_applets.inp']);
    reversibility_test = InputFile([path_inp_files,'/reversibility_test.inp']);

    %read default inputfiles
    %All variables get a default value from the blueprints
    gorilla.read();
    tetra_grid.read();
    gorilla_applets.read();
    reversibility_test.read();


    %Set input variables for GORILLA
    %-----------------------------------------------------------------------------------------------------------------------

    %Input file gorilla.inp
    %All necessary variables for current calculation

        %Change in the electrostatic potential within the plasma volume in Gaussian units
            gorilla.GORILLANML.eps_Phi = eps_Phi;

        %Coordinate system
            %1 ... (R,phi,Z) cylindrical coordinate system
            %2 ... (s,theta,phi) symmetry flux coordinate system
            gorilla.GORILLANML.coord_system = 2;

        %particle species
            %1 ... electron, 2 ... deuterium ion, 3 ... alpha particle
            gorilla.GORILLANML.ispecies = ispecies;

        %Switch for initial periodic coordinate particle re-location (modulo operation)
            % true ... Particles are re-located at initialization in the case of a periodic coordinate, if they are outside the computation domain.
            % false ... Particles are not re-located at initialization (This might lead to error if particles are outside the computation domain)
            gorilla.GORILLANML.boole_periodic_relocation = true;

        %1 ... numerical RK pusher, 2 ... polynomial pusher
            gorilla.GORILLANML.ipusher = 2;

        %true ... Face guessing algorithm is used, false ... NOT used
            gorilla.GORILLANML.boole_guess = boole_guess;

        %Polynomial order for orbit pusher (from 2 to 4)
            gorilla.GORILLANML.poly_order = poly_order;

        %Time tracing options
        % 1 ... Time tracing in 0th order approximation: dt / dtau = const. in each tetrahedral cell
        % 2 ... (Correct) Hamiltonian time tracing in adequate polynomial order
        %       Dependencies:  ipusher must be 2 (Polynomial pusher)
            gorilla.GORILLANML.i_time_tracing_option = i_time_tracing_option;
            
        % Switches to calculate optional quantities
            gorilla.GORILLANML.boole_time_Hamiltonian = true;
            gorilla.GORILLANML.boole_gyrophase = true;

        % Switch for adaptive time step scheme
            % false ... no adaptive scheme
            % true ... adaptive scheme to ensure energy conservation up to specified fluctuation
            gorilla.GORILLANML.boole_adaptive_time_steps = false;


    %Input file tetra_grid.inp
    %All necessary variables for current calculation

        %Grid Size
            %Rectangular: nR, Field-aligned: ns
            tetra_grid.TETRA_GRID_NML.n1 = n1;
            %Rectangular: nphi, Field-aligned: nphi
            tetra_grid.TETRA_GRID_NML.n2 = n2;
            %Rectangular: nZ, Field-aligned: ntheta
            tetra_grid.TETRA_GRID_NML.n3 = n3;

        %Grid kind
            %1 ... rectangular grid for axisymmetric EFIT data
            %2 ... field-aligned grid for axisymmetric EFIT data
            %3 ... field-aligned grid for non-axisymmetric VMEC
            %4 ... SOLEDGE3X_EIRENE grid
            tetra_grid.TETRA_GRID_NML.grid_kind = grid_kind;

        %Switch for selecting number of field periods automatically or manually
            %.true. ... number of field periods is selected automatically (Tokamak = 1, Stellarator depending on VMEC equilibrium)
            %.false. ... number of field periods is selected manually (see below)
            tetra_grid.TETRA_GRID_NML.boole_n_field_periods = true;
            
        %Symmetry Flux Coordinates Annulus (Minimal value for flux coordinates s)
            tetra_grid.TETRA_GRID_NML.sfc_s_min = sfc_s_min;


    %Input file gorilla_applets.inp
    %All necessary variables for current calculation

        % Choice of application
            % 1 ... Pre-computation of fluxtube volume for GORILLA
            % 2 ... Computation of mono-energetic transport coefficient for given collisionality and Mach number
            % 3 ... Scan over collisionality: Computation of mono-energetic transport coefficent for constant Mach number
            % 4 ... Computation of numerical diffusion coefficient (No collisions)
            % 5 ... Computation of alpha particle lifetime
            % 6 ... Demonstrate direct VMEC integrator
            % 7 ... Computation of first Poincare invariance
            % 8 ... Perform a reversibility test (demonstrate negative time)
            % 9 ... Computation of particle dwell times in tetrahedron cells
            % 10 .. Perform a bolzmann test 
            gorilla_applets.GORILLA_APPLETS_NML.i_option = 8;
            
            
    %Input file reversibility_test.inp
    %All necessary variables for current calculation

        %Total time for tracing particles
        % HYDRA ions (3keV);        t = 1.3d-4
        % TOk ions (3keV);          t = 3.25d-5
        % TOK electrons (1E-2eV)    t = 2.1d-4
            reversibility_test.REVERSIBILITY_TEST_NML.t_total = t_total;

        %Average "minor radius" of orbits at start
        %Cylindrical: r, Symmetry Flux: s
            reversibility_test.REVERSIBILITY_TEST_NML.minor_0 = 0.5d0;
            
        %Bandwith of "minor radius" of orbits at start
        %Cylindrical: dr, Symmetry Flux: ds
            reversibility_test.REVERSIBILITY_TEST_NML.dminor_0 = dminor_0;

        %Major Radius R of starting contour center
        %ONLY NEEDED FOR coord_system = 1 (cylindrical coordinates)
            reversibility_test.REVERSIBILITY_TEST_NML.R_center = 2.4d0;

        %Height Z of starting contour center
        %ONLY NEEDED FOR coord_system = 1 (cylindrical coordinates)
            reversibility_test.REVERSIBILITY_TEST_NML.Z_center = 0.0d0;
            
        %Fraction of the contour to be traced
            reversibility_test.REVERSIBILITY_TEST_NML.contour_fraction = contour_fraction;

        %Starting $\varphi$ value of orbits at start
            reversibility_test.REVERSIBILITY_TEST_NML.phi_0 = phi_0;

        %Average starting pitch parameter of particles
            reversibility_test.REVERSIBILITY_TEST_NML.pitchpar_0 = pitchpar_0;
            
        %Bandwith of starting pitch parameter of particles
            reversibility_test.REVERSIBILITY_TEST_NML.dpitchpar = dpitchpar;

        %Average kinetic energy of particles at start
            reversibility_test.REVERSIBILITY_TEST_NML.energy_eV_0 = energy_eV_0;

        %RELATIVE Bandwith of energy of particles at start
            reversibility_test.REVERSIBILITY_TEST_NML.relative_bandwith = relative_bandwith;

        %Number of time measurements
            reversibility_test.REVERSIBILITY_TEST_NML.n_steps = n_steps;

        %Number of orbits
            reversibility_test.REVERSIBILITY_TEST_NML.n_orbits = n_orbits;

        %Number of time snapeshots
            reversibility_test.REVERSIBILITY_TEST_NML.n_snapshots = n_snapshots;

        %Filename for reversibility test (forward part)
            reversibility_test.REVERSIBILITY_TEST_NML.filename_reversibility_test = filename_reversibility_test;
            
        %Filename for reversibility test (backward part)
            reversibility_test.REVERSIBILITY_TEST_NML.filename_reversibility_test_back = filename_reversibility_test_back;

        %Number of time measurements
            reversibility_test.REVERSIBILITY_TEST_NML.boole_diag_reversibility_test = boole_diag_reversibility_test;


    %Run GORILLA
    %-----------------------------------------------------------------------------------------------------------------------

    %write Input files for GORILLA
    gorilla.write([path_RUN,'/gorilla.inp']);
    tetra_grid.write([path_RUN,'/tetra_grid.inp']);
    reversibility_test.write([path_RUN,'/reversibility_test.inp']);
    gorilla_applets.write([path_RUN,'/gorilla_applets.inp']);

    %Create softlinks for used files
    if isfile('../../../gorilla_applets_main.x')
        ! ln -s ../../../gorilla_applets_main.x .
    elseif isfile('../../../BUILD/SRC/gorilla_applets_main.x')
        ! ln -s ../../../BUILD/SRC/gorilla_applets_main.x .
    else
        disp('GORILLA_APPLETS not built, exiting the MatLab script')
        return
    end
    ! ln -s ../../../../GORILLA/MHD_EQUILIBRIA .
    ! ln -s ../../../../GORILLA/INPUT/field_divB0.inp .
    ! ln -s ../../../../GORILLA/INPUT/preload_for_SYNCH.inp .

    %Run without energy control
    if (boole_recalculate)
        ! ./gorilla_applets_main.x
    end
    
    cd(path_script)
    
end

%Create plots of generated data
%-----------------------------------------------------------------------------------------------------------------------

%Load results
if (boole_show_contour)
    contour_data = load([path_RUN,'/',filename_reversibility_test]);
end
if (boole_show_contour_back)
    contour_back_data = load([path_RUN,'/',filename_reversibility_test_back]);
end

%Show contour evolution
if (boole_show_contour)
    xlabel_txt = '$k_{orbit}$';
    
    ylimits = nan*ones(2,number_of_quantities);
    [~,contour_data] = extract_snapshot(contour_data);
    ylimits(1,:) = min(contour_data(:,2:number_of_quantities+1),[],1);
    ylimits(2,:) = max(contour_data(:,2:number_of_quantities+1),[],1);
    
    figure('Renderer', 'painters', 'Position', [24 37 2500 1300])
    t = tiledlayout(number_of_quantities,n_snapshots +1);
    
    for j = 1:(n_snapshots+1)
        plot_data = extract_snapshot(contour_data,j);
        for l = 1:number_of_quantities
            jl_entry = (l-1)*(n_snapshots+1) + j;
            contour_labels = [1:size(plot_data,1)];
            %subplot(number_of_quantities,n_snapshots + 1,jl_entry)
            nexttile(jl_entry)
            plot(contour_labels(1:n_skip_contour:end),plot_data(1:n_skip_contour:end,l+1),MarkerContour,'Color',forwardColor,'MarkerSize',MarkerSizeContour)
            ylim(ylimits(:,l))
            if (j == 1 || l == 1 || l == number_of_quantities)
                ax = gca;
                ax.FontSize = FontSize*axis_factor;
                if (l == 1)
                    tit = title(['$\tau = $ ',num2str((j-1)/n_snapshots)],'Interpreter','latex');
                    tit.FontSize = FontSize;
                end
                if (l == number_of_quantities)
                    xlab = xlabel(xlabel_txt,'Interpreter','latex');
                    xlab.FontSize = FontSize;
                else
                    set(gca,'XTick',[])
                end
                if (j == 1)
                    ylabel_txt = ylabel_quantities{l};
                    ylab = ylabel(ylabel_txt,'Interpreter','latex');
                    ylab.FontSize = FontSize;
                else
                    set(gca,'YTick',[])
                end
            else
                set(gca,'XTick',[])
                set(gca,'YTick',[])
            end
        end
    end
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
end

%Show contour evolution backwards
if (boole_show_contour_back)
    xlabel_txt = '$k_{orbit}$';
    
%     ylimits = nan*ones(2,number_of_quantities);
%     [~,contour_back_data] = extract_snapshot(contour_back_data);
%     ylimits(1,:) = min(contour_back_data(:,2:number_of_quantities+1),[],1);
%     ylimits(2,:) = max(contour_back_data(:,2:number_of_quantities+1),[],1);
    
    figure('Renderer', 'painters', 'Position', [24 37 2500 1300])
    t = tiledlayout(number_of_quantities,n_snapshots +1);
    
    for j = 1:(n_snapshots+1)
        plot_data = extract_snapshot(contour_back_data,j);
        for l = 1:number_of_quantities
            jl_entry = (l-1)*(n_snapshots+1) + j;
            contour_labels = [1:size(plot_data,1)];
            %subplot(number_of_quantities,n_snapshots + 1,jl_entry)
            nexttile(jl_entry)
            plot(contour_labels(1:n_skip_contour:end),plot_data(1:n_skip_contour:end,l+1),MarkerContour,'Color',forwardColor,'MarkerSize',MarkerSizeContour)
            ylim(ylimits(:,l))
            if (j == 1 || l == 1 || l == number_of_quantities)
                ax = gca;
                ax.FontSize = FontSize*axis_factor;
                if (l == 1)
                    tit = title(['$\tau = $ ',num2str((j-1)/n_snapshots)],'Interpreter','latex');
                    tit.FontSize = FontSize;
                end
                if (l == number_of_quantities)
                    xlab = xlabel(xlabel_txt,'Interpreter','latex');
                    xlab.FontSize = FontSize;
                else
                    set(gca,'XTick',[])
                end
                if (j == 1)
                    ylabel_txt = ylabel_quantities{l};
                    ylab = ylabel(ylabel_txt,'Interpreter','latex');
                    ylab.FontSize = FontSize;
                else
                    set(gca,'YTick',[])
                end
            else
                set(gca,'XTick',[])
                set(gca,'YTick',[])
            end
        end
    end
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HELPFUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [snapshot,array] = extract_snapshot(array,order)

    I = isnan(array(:,1));
    if ~any(I)
        I = array(:,1) == -1;
        array(I,:) = nan;
        if isempty(I)
            error('Error in extract_instance: Array does not contain any known seperators (-1 or nan)!');
        end
    end
    index_set = find(I) + 1;

    if (nargin < 2 || isempty(order))
        order = length(index_set);
    end
    if (order < 1 || order > length(index_set))
        error(['Error in extract_instance: Requested target value is out of range [', num2str(1)',',',num2str(length(index_set)),']!']);
    end
    if (order == 1)
        index_snapshot_start = 1;
    else
        index_snapshot_start = index_set(order-1);
    end
    index_snapshot_end = index_set(order) - 2;
    
    snapshot_indices = [index_snapshot_start,index_snapshot_end];
    snapshot = array(snapshot_indices(1):snapshot_indices(2),:);
end
%%
% % Load Poincaré invariant over time steps
% poincare_invariant_data = load(['./results/poincare_invariance_gorilla_n_1E7_e_var_jperp_var_dtdtau_ham.dat']);
% 
% 
% 
% f1 = figure;
% f1.Units = 'normalized';
% f1.Position = [0.0797 0.1657 0.6271 0.6565];  
%     
% %
% p_1 = plot(poincare_invariant_data(:,1),poincare_invariant_data(:,2),'-');
% p_1.DisplayName = '$a_\perp = 0.1$, $a_H = 0.1$: Correct (Hamiltonian) gyro-phase $\phi$ \& time dynamics $\mathrm{d} t / \mathrm{d} \tau$';
% 
% %
% l_1 = legend;
% set([p_1],'LineWidth',3);
% set(gca,'FontSize',18);
% l_1.FontSize = 22;
% l_1.Location = 'northwest';
% l_1.Interpreter = 'latex';
% %
% ylabel('$J$','interp','latex')
% xlabel('$\tau_\mathrm{step}$','interp','latex')
% %
% box on
% xlim([0,200]);
% 
% 
% 
% 
% 
% f2 = figure;
% f2.Units = 'normalized';
% f2.Position = [0.0797 0.1657 0.6271 0.6565];   
% %
% p_2 = semilogy(normalize_tau(poincare_invariant_data(:,1)),normalize_pi(poincare_invariant_data));
% p_2.DisplayName = p_1.DisplayName;
% p_2.Color = p_1.Color;
% p_2.Marker = p_1.Marker;
% p_2.LineStyle = p_1.LineStyle;
% 
% %
% l_2 = legend;
% set([p_2],'LineWidth',3);
% set(gca,'FontSize',18);
% l_2.Location = 'southwest';
% l_2.FontSize = 22;
% l_2.Interpreter = 'latex';
% %
% ylabel('$\Delta J ~/~ J (\tau=0)$','interp','Latex')
% xlabel('$\tau~/~\tilde{\tau}$','interp','latex')
% %
% box on
% grid on
% xlim([0,1]);
% ylim([1E-20,1E-2]);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% function [inv_out] = normalize_pi(inv_in)
%     inv_out(:) = abs((inv_in(:,2)-inv_in(1,2))/inv_in(1,2));
% end
% 
% function [tau_out] = normalize_tau(tau_in)
%     tau_out(:) = tau_in(:)/tau_in(end);
% end
