%#######################################################################################################################
% description:
%-----------------------------------------------------------------------------------------------------------------------
% * Compute collisionless guiding-center orbit with GORILLA for a contour of Deuterium particles with various starting conditions.
% * Use a field-aligned grid for several MHD equilibria (ASDEX, VMEC, WEST).
% * Perform a reversibility test by integrating back from the end postion of the orbit.
% * Compare the forward/backward evolution of the contour via snapshots.
%
%#######################################################################################################################
% used functions:
%-----------------------------------------------------------------------------------------------------------------------
% InputFile
% NameList
% read_in
%#######################################################################################################################
% authors: Georg Graßler
% created: 10.10.2022

%For Construction of Demonstration
%-----------------------------------------------------------------------------------------------------------------------
%control of script
close all
boole_recalculate = false;
ipusher = 2;
poly_order = 2;
grid_kind = 4;
coord_system = 2;
[n1,n2,n3] = deal(70,70,70);
[s_0,ds_0,theta_0,dtheta_0,R_0,dR_0,Z_0,dZ_0] = deal(-1);
[g_file_filename,convex_wall_filename] = deal('x');
sfc_s_min = 1.d-3;
eps_Phi = -1.d-7;
%eps_Phi = 0;
ispecies = 2;
boole_guess = true;
i_time_tracing_option = 1;
contour_fraction = 1;
n_steps = 20;
n_orbits = 99;
n_snapshots = 2;
energy_eV_0 = 3d3;
boole_apply_noise = false;
seed_option = 2;
boole_diag_reversibility_test = false;

boole_long_test = false;

switch(grid_kind)
    case(2)
        t_total = 3.25d-5*6/7;
        %t_total = 2.0d-4;
        %t_total = 7d-5;
        if (boole_long_test)
            n_steps = 1;
            n_snapshots = 1;
            t_total = t_total*1000;
        end
        s_0 = 0.5d0;
        ds_0 = 0.25d0;
        theta_0 = 1.5d0;
        dtheta_0 = 1.0d0;
        phi_0 = 0.01;
        pitchpar_0 = 0.46d0;
        dpitchpar = 0.3d0;
        relative_bandwith = 0.1d0;
        noise_amplitude = 0.02d0;
        g_file_filename = 'MHD_EQUILIBRIA/g_file_for_test';
        convex_wall_filename = 'MHD_EQUILIBRIA/convex_wall_for_test.dat';
        ylabel_quantities = {'$t$','$s$','$\vartheta$','$\varphi$','$\lambda$','$\phi$','$E_{kin}$'};
    case(3)
        t_total = 1.3d-4*3/4;
        %t_total = 2d-4;
        coord_system = 2;
        eps_Phi = -5.d-6;
        s_0 = 0.5d0;
        ds_0 = 0.4d0;
        theta_0 = 5.0d0;
        %theta_0 = 0.0d0;
        dtheta_0 = 1.0d0;
        phi_0 = 1.0d0;
        pitchpar_0 = 0.4d0;
        dpitchpar = -0.3d0;
        relative_bandwith = 0.9d0;
        noise_amplitude = 0.01d0;
        ylabel_quantities = {'$t$','$s$','$\vartheta$','$\varphi$','$\lambda$','$\phi$','$E_{kin}$'};
    case(4)
        t_total = 4.2d-5;
        t_total = 8.2d-5;
        coord_system = 1;
        eps_Phi = +1.d-6;
        dR_0 = 1.0d0;
        dZ_0 = 1.0d0;
%         R_0 = 231.0d0;
%         R_0 = 245.0d0;
        R_0 = 269.0d0;
%         Z_0 = -53.0d0;
%         Z_0 = -51.0d0;
        Z_0 = -40.0d0;
%         dR_0 = 1.5d0;
%         dZ_0 = 1.5d0;
        phi_0 = 5.0d0;
        pitchpar_0 = -0.3d0;
        dpitchpar = -0.2d0;
        relative_bandwith = 0.1d0;
        noise_amplitude = 0.01d0;
        n3 = 30;
        g_file_filename = 'MHD_EQUILIBRIA/g_file_for_test_WEST';
        convex_wall_filename = 'MHD_EQUILIBRIA/convex_wall_for_test_WEST.dat';
        ylabel_quantities = {'$t$','$R$','$\varphi$','$Z$','$\lambda$','$\phi$','$E_{kin}$'};
end

boole_show_contour = false;
boole_show_contour_back = true;
boole_show_contour_diff = true;
boole_show_WEST_contour = true;
number_of_quantities = 7;
has_central_line = {'$t$','$\lambda$'};
has_normalized_limits = {'$\lambda$'};
n_skip_contour = 1;

forwardColor = [215,25,28]/256;
backwardColor = [44,123,182]/256;
DiffColor = [0,0,0]/256;
CentralLineColor = [100,100,100]/256;

MarkerContour = '.';
MarkerContourBack = 'o';
MarkerDiff = '.';
MarkerSizeContour = 8;
MarkerSizeContourBack = 5;
MarkerSizeDiff = 8;
CentralLineWidth = 2;
FontSize = 25;
FontSizeTitle = 20;
axis_factor = 3/5;

postprocessing_functions_path = 'postprocessing_functions';
addpath(postprocessing_functions_path);

%filenames for output
filename_reversibility_test = 'reversibility_test.dat';
filename_reversibility_test_back = 'reversibility_test_back.dat';

  
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

if (boole_recalculate)
    
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
            gorilla.GORILLANML.coord_system = coord_system;

        %particle species
            %1 ... electron, 2 ... deuterium ion, 3 ... alpha particle
            gorilla.GORILLANML.ispecies = ispecies;

        %Switch for initial periodic coordinate particle re-location (modulo operation)
            % true ... Particles are re-located at initialization in the case of a periodic coordinate, if they are outside the computation domain.
            % false ... Particles are not re-located at initialization (This might lead to error if particles are outside the computation domain)
            gorilla.GORILLANML.boole_periodic_relocation = true;

        %1 ... numerical RK pusher, 2 ... polynomial pusher
            gorilla.GORILLANML.ipusher = ipusher;

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
            
        %MHD equilibrium filename
            tetra_grid.TETRA_GRID_NML.g_file_filename = g_file_filename;
            tetra_grid.TETRA_GRID_NML.convex_wall_filename = convex_wall_filename;
            tetra_grid.TETRA_GRID_NML.netcdf_filename = 'MHD_EQUILIBRIA/netcdf_file_for_test.nc';

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

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %ONLY FOR coord_system = 1 (cylindrical coordinates)
        %Fluxlabel s of starting contour center
        reversibility_test.REVERSIBILITY_TEST_NML.s_0 = s_0;

        %Bandwith of flux label s of starting contour
        reversibility_test.REVERSIBILITY_TEST_NML.ds_0 = ds_0;

        %Theta of starting contour center
        reversibility_test.REVERSIBILITY_TEST_NML.theta_0 = theta_0;

        %Bandwith of theta of starting contour
        reversibility_test.REVERSIBILITY_TEST_NML.dtheta_0 = dtheta_0;
            
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %ONLY FOR coord_system = 1 (cylindrical coordinates)
        %Major Radius R of starting contour center
            reversibility_test.REVERSIBILITY_TEST_NML.R_0 = R_0;

        %Bandwitch of radius R of starting contour
            reversibility_test.REVERSIBILITY_TEST_NML.dR_0 = dR_0;

        %Height Z of starting contour center
            reversibility_test.REVERSIBILITY_TEST_NML.Z_0 = Z_0;

        %Bandwitch of height Z of starting contour
             reversibility_test.REVERSIBILITY_TEST_NML.dZ_0 = dZ_0;
             
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Starting $\varphi$ value of orbits at start
            reversibility_test.REVERSIBILITY_TEST_NML.phi_0 = phi_0;
            
        %Fraction of the contour to be traced
            reversibility_test.REVERSIBILITY_TEST_NML.contour_fraction = contour_fraction;

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
            
        %Switch to apply noise on coordinates at reversal point
            reversibility_test.REVERSIBILITY_TEST_NML.boole_apply_noise = boole_apply_noise;

        %Realtive amplitude of noise on coordinates at reversal point
            reversibility_test.REVERSIBILITY_TEST_NML.noise_amplitude = noise_amplitude;

        %Option for computation of above random noise
        % 1 ... compute seed with a random number
        % 2 ... load seed from file 'seed.inp'
            reversibility_test.REVERSIBILITY_TEST_NML.seed_option = seed_option;

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
if (boole_show_contour || boole_show_contour_back || boole_show_WEST_contour)
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

a = extract_snapshot(contour_data,2);

%Show contour evolution backwards
if (boole_show_contour_back)
    xlabel_txt = '$k_{orbit}$';
    
    ylimits = nan*ones(2,number_of_quantities);
    ylimits_back = nan*ones(2,number_of_quantities);
    [~,contour_data] = extract_snapshot(contour_data);
    [~,contour_back_data] = extract_snapshot(contour_back_data);
    ylimits(1,:) = min(contour_data(:,2:number_of_quantities+1),[],1);
    ylimits(2,:) = max(contour_data(:,2:number_of_quantities+1),[],1);
    ylimits_back(1,:) = min(contour_back_data(:,2:number_of_quantities+1),[],1);
    ylimits_back(2,:) = max(contour_back_data(:,2:number_of_quantities+1),[],1);
    ylimits(1,:) = min([ylimits(1,:);ylimits_back(1,:)],[],1);
    ylimits(2,:) = max([ylimits(2,:);ylimits_back(2,:)],[],1);
    
    figure('Renderer', 'painters', 'Position', [24 37 2500 1300])
    t = tiledlayout(number_of_quantities,n_snapshots +1);
    
    for j = 1:(n_snapshots+1)
        plot_data = extract_snapshot(contour_data,j);
        plot_data_back = extract_snapshot(contour_back_data,j);
        for l = 1:number_of_quantities
            jl_entry = (l-1)*(n_snapshots+1) + j;
            ylabel_txt = ylabel_quantities{l};
            contour_labels = [1:size(plot_data_back,1)];
            %subplot(number_of_quantities,n_snapshots + 1,jl_entry)
            nexttile(jl_entry)
            if(ismember(ylabel_txt,has_central_line))
                plot(contour_labels,zeros(size(contour_labels)),'-','Color',CentralLineColor,'LineWidth',CentralLineWidth)
            end
            hold on
            plot(contour_labels(1:n_skip_contour:end),plot_data(1:n_skip_contour:end,l+1),MarkerContour,'Color',forwardColor,'MarkerSize',MarkerSizeContour,'MarkerEdgeColor',forwardColor)
            plot(contour_labels(1:n_skip_contour:end),plot_data_back(1:n_skip_contour:end,l+1),MarkerContourBack,'Color',backwardColor,'MarkerSize',MarkerSizeContourBack,'MarkerEdgeColor',backwardColor)
            hold off
            if(ismember(ylabel_txt,has_normalized_limits))
                ylim([-1,1])
            else
                ylim(ylimits(:,l))
            end
            xlim([contour_labels(1),contour_labels(end)])
            if (j == 1 || l == 1 || l == number_of_quantities)
                ax = gca;
                ax.FontSize = FontSize*axis_factor;
                if (l == 1)
                    %tit = title(['$\tau = $ ',num2str((j-1)/n_snapshots)],'Interpreter','tex');
                    forward_title = ['\tau_{+} = ',num2str((j-1)/n_snapshots)];
                    backward_title = ['\tau_{-} = ',num2str((n_snapshots + 1 - j)/n_snapshots)];
                    tit = title(sprintf('\\color[rgb]{%f, %f, %f}%s \\color{black}/ \\color[rgb]{%f, %f, %f}%s', forwardColor, forward_title,backwardColor,backward_title),'Interpreter','tex');
                    tit.FontSize = FontSizeTitle;
                end
                if (l == number_of_quantities)
                    xlab = xlabel(xlabel_txt,'Interpreter','latex');
                    xlab.FontSize = FontSize;
                else
                    set(gca,'Xticklabel',[])
                end
                if (j == 1)
                    ylab = ylabel(ylabel_txt,'Interpreter','latex');
                    ylab.FontSize = FontSize;
                else
                    set(gca,'Yticklabel',[])
                end
            else
                set(gca,'Xticklabel',[])
                set(gca,'Yticklabel',[])
            end
        end
    end
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
end

%Show contour evolution backwards
if (boole_show_contour_diff)
    xlabel_txt = '$k_{orbit}$';
    
    ylimits = nan*ones(2,number_of_quantities);
    ylimits_back = nan*ones(2,number_of_quantities);
    ylimits_diff = nan*ones(2,number_of_quantities);
    [~,contour_data] = extract_snapshot(contour_data);
    [~,contour_back_data] = extract_snapshot(contour_back_data);
    ylimits(1,:) = min(contour_data(:,2:number_of_quantities+1),[],1);
    ylimits(2,:) = max(contour_data(:,2:number_of_quantities+1),[],1);
    ylimits_back(1,:) = min(contour_back_data(:,2:number_of_quantities+1),[],1);
    ylimits_back(2,:) = max(contour_back_data(:,2:number_of_quantities+1),[],1);
    ylimits(1,:) = min([ylimits(1,:);ylimits_back(1,:)],[],1);
    ylimits(2,:) = max([ylimits(2,:);ylimits_back(2,:)],[],1);
    
%     contour_diff_data = contour_data;
%     contour_diff_data(:,2:end) = (contour_data(:,2:end)-contour_back_data(:,2:end))./max(abs(contour_data(:,2:end)));
%     ylimits_diff(1,:) = min(contour_diff_data(:,2:number_of_quantities+1),[],1);
%     ylimits_diff(2,:) = max(contour_diff_data(:,2:number_of_quantities+1),[],1);
    [contour_diff_data,ylimits_diff] = get_diff(contour_data,contour_back_data,number_of_quantities);

    figure('Renderer', 'painters', 'Position', [24 37 2500 1300])
    t = tiledlayout(number_of_quantities,n_snapshots + 1);
    
    for j = [1:(n_snapshots+1)]
        plot_data = extract_snapshot(contour_data,j);
        plot_data_back = extract_snapshot(contour_back_data,j);
        plot_data_diff = extract_snapshot(contour_diff_data,j);
        for l = 1:number_of_quantities
            jl_entry = (l-1)*(n_snapshots + 1) + j;
            ylabel_txt = ylabel_quantities{l};
            contour_labels = [1:size(plot_data_back,1)];
            %subplot(number_of_quantities,n_snapshots + 1,jl_entry)
            nexttile(jl_entry)
            hold on
            if(ismember(ylabel_txt,['$t$']))
                plot(contour_labels,zeros(size(contour_labels)),'-','Color',CentralLineColor,'LineWidth',CentralLineWidth)
                hold on
                plot(contour_labels(1:n_skip_contour:end),plot_data(1:n_skip_contour:end,l+1),MarkerContour,'Color',forwardColor,'MarkerSize',MarkerSizeContour,'MarkerEdgeColor',forwardColor)
                plot(contour_labels(1:n_skip_contour:end),plot_data_back(1:n_skip_contour:end,l+1),MarkerContourBack,'Color',backwardColor,'MarkerSize',MarkerSizeContourBack,'MarkerEdgeColor',backwardColor)
                hold off
                ylim(ylimits(:,l))
            else
                hold on
                plot(contour_labels(1:n_skip_contour:end),plot_data_diff(1:n_skip_contour:end,l+1),MarkerDiff,'Color',DiffColor,'MarkerSize',MarkerSizeDiff,'MarkerEdgeColor',DiffColor)
                hold off
                ylim(ylimits_diff(:,l))
                grid on
            end
            xlim([contour_labels(1),contour_labels(end)])
            if (j == 1 || l == 1 || l == number_of_quantities)
                ax = gca;
                ax.FontSize = FontSize*axis_factor;
                if (l == 1)
                    forward_title = ['\tau_{+} = ',num2str((j-1)/n_snapshots)];
                    backward_title = ['\tau_{-} = ',num2str((n_snapshots + 1 - j)/n_snapshots)];
                    tit = title(sprintf('\\color[rgb]{%f, %f, %f}%s \\color{black}/ \\color[rgb]{%f, %f, %f}%s', forwardColor, forward_title,backwardColor,backward_title),'Interpreter','tex');
                    tit.FontSize = FontSizeTitle;
                end
                if (l == number_of_quantities)
                    xlab = xlabel(xlabel_txt,'Interpreter','latex');
                    xlab.FontSize = FontSize;
                else
                    set(gca,'Xticklabel',[])
                end
                if (j == 1)
                    ylab = ylabel(ylabel_txt,'Interpreter','latex');
                    ylab.FontSize = FontSize;
                else
                    set(gca,'Yticklabel',[])
                end
            else
                set(gca,'Xticklabel',[])
                set(gca,'Yticklabel',[])
            end
        end
    end
    
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
end

if (grid_kind == 4 && boole_show_WEST_contour)
    %Color scheme
    grid_color = [204,204,204]/256;
    grid_thickness = 0.5;
    WESTColor = [4,90,141]/256;
    MarkerWEST = '.';
    MarkerSizeWEST = 13;
    xlabel_txt = '$R$ [cm]';
    ylabel_txt = '$Z$ [cm]';

    %Read in Soledge3x-EIRENE mesh data and prepare 2D grid
    fileID = fopen([path_main_core,'/MHD_EQUILIBRIA/MESH_SOLEDGE3X_EIRENE/knots_for_test.dat']);
    coordinates = textscan(fileID,'%f %f','HeaderLines',1);
    fclose(fileID);
    coordinates = cell2mat(coordinates);
    n_vertex = length(coordinates);

    fileID = fopen([path_main_core,'/MHD_EQUILIBRIA/MESH_SOLEDGE3X_EIRENE/triangles_for_test.dat']);
    triangles = textscan(fileID,'%f %f %f','HeaderLines',1);
    fclose(fileID);
    triangles = cell2mat(triangles);
    n_triangles = length(triangles);

    grid = nan*ones(n_triangles*5,2);
    for t = [1:n_triangles]
        for k = [1:3]
            grid((t-1)*5 + k,:) = coordinates(triangles(t,k),:);
        end
        grid((t-1)*5 + 4,:) = coordinates(triangles(t,1),:);
        grid((t-1)*5 + 5,:) = nan;
    end

    figure('Renderer', 'painters', 'Position', [100 100 1300 1000])
    grid_plot = plot(grid(:,1),grid(:,2),'Color',grid_color);
    grid_plot.LineWidth = grid_thickness;
    grid_plot.DisplayName = 'SOLEDGE3X-EIRENE mesh';
    contour_handles = cell(1,n_snapshots+1);
    hold on
    for j = [1:1:n_snapshots+1]
        plot_data = extract_snapshot(contour_data,j);
        %plot(plot_data(:,3),plot_data(:,5),MarkerWEST,'Color',WESTColor,'MarkerSize',MarkerSizeWEST)
        contour_handles{j} = plot(plot_data(:,3),plot_data(:,5),MarkerWEST,'MarkerSize',MarkerSizeWEST);
        contour_handles{j}.DisplayName = [num2str(j),'. Snapshot'];
    end
    hold off
    len = legend([grid_plot,contour_handles{:}],'Interpreter','latex');
    xlab = xlabel(xlabel_txt,'Interpreter','latex');
    ylab = ylabel(ylabel_txt,'Interpreter','latex');
    ax = gca;
    ax.FontSize = FontSize*axis_factor;
    xlab.FontSize = FontSize;
    ylab.FontSize = FontSize;
    len.FontSize = FontSize;
    
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

function [array_diff,ylimits_diff] = get_diff(array,array_back,number_of_quantities)
    array_diff = array;
%     distance = max(abs(extract_snapshot(array,1)-extract_snapshot(array)));
%     distance = distance(2:end);
    distance = max(abs(array(:,2:end)));
    diff = abs(array(:,2:end)-array_back(:,2:end));
    array_diff(:,2:end) = diff./distance;
    ylimits_diff(1,:) = min(array_diff(:,2:number_of_quantities+1),[],1);
    ylimits_diff(2,:) = max(array_diff(:,2:number_of_quantities+1),[],1);
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
