% Precomputation of alpha lifetime computation


% Save MATLAB path
matlab_path = pwd;

% local path of libneo
% libneo_path = '~/GITHUB/libneo';

% local path of GORILLA
gorilla_path = [matlab_path,'/../'];  % '~/GITHUB/gorilla/';

% add path of locale installation of InputFile class
addpath([gorilla_path, 'MATLAB/functions']);

 
% initialize input file class
gorilla = InputFile([gorilla_path,'inp_blueprints/gorilla.inp']);
gorilla_applets = InputFile([gorilla_path,'inp_blueprints/gorilla_applets.inp']);
tetra_grid = InputFile([gorilla_path,'inp_blueprints/tetra_grid.inp']);
alpha_lifetime = InputFile([gorilla_path,'inp_blueprints/alpha_lifetime.inp']);
  
% read file
gorilla.read();
gorilla_applets.read();
tetra_grid.read();
alpha_lifetime.read();

% General settings for Precomputation of Alpha Lifetime Runs

% Name of path with Fluxtubevolumes
folder_name_fluxtv = 'RUN/fluxtubevolumes';


%% General settings for computation of Alpha Lifetime Runs

% Set computation mode to FluxTV computation
gorilla_applets.GORILLA_APPLETS_NML.i_option = 5;

% Choose alphas as species
gorilla.GORILLANML.ispecies = 3;

% Pusher type
gorilla.GORILLANML.ipusher = 2;

% Poly order
gorilla.GORILLANML.poly_order = 2;

% Grid kind (Stellarator)
grid_kind = 3;
tetra_grid.TETRA_GRID_NML.grid_kind = grid_kind;

% Coordinate system (Symmetry Flux Coordinates)
coord_system = 2;
gorilla.GORILLANML.coord_system = coord_system;

% Alpha lifetime settings
%Set integrator type
% 1... GORILLA
% 2 ... DIRECT
alpha_lifetime.alpha_lifetimenml.i_integrator_type = 1;
alpha_lifetime.alpha_lifetimenml.time_step = 0.01;
alpha_lifetime.alpha_lifetimenml.n_particles = 10000;


%% Computation of flux tube volumes

% Create RUN folder for Alpha Lifetime Runs
name_run_folder = 'RUN/20220816_alpha_lifetime_1.d-2s_Ns_100_Nang_scan';
mkdir(gorilla_path,name_run_folder);

% Create directory for results
mkdir([gorilla_path,name_run_folder,'/'],'results');


%% Loop over


%grid_vec = [28,39,53,73,101,139,191,264,363,500];
grid_vec = [28,32];

[NR] = deal(100);

%grid_vec = [70,100,200];

%[NPHI,NTHETA] = deal(100);


for i = 1:numel(grid_vec)

    % Create directory for run
    folder_name = ['folder_Poly4_',num2str(grid_vec(i))];
    %folder_name = 'folder_odeint_1dm8';

    mkdir([gorilla_path,name_run_folder,'/'],folder_name);

    % Create softlinks
    status = system(['ln -s ',gorilla_path,'gorilla_applets_main.x ',gorilla_path,name_run_folder,'/',folder_name,'/.'],'-echo')
    status = system(['ln -s ',gorilla_path,'inp_blueprints/seed.inp ',gorilla_path,name_run_folder,'/',folder_name,'/.'],'-echo')
    status = system(['ln -s ',gorilla_path,'DATA ',gorilla_path,name_run_folder,'/',folder_name,'/.'],'-echo')
    status = system(['ln -s ',gorilla_path,'DATA/VMEC/wout_23_1900_fix_bdry.nc ',gorilla_path,name_run_folder,'/',folder_name,'/netcdf_file_for_test.nc'],'-echo')

    % Set grid size
    [NPHI,NTHETA] = deal(grid_vec(i));

    % Define filename and path
    file_name_fluxtv = ['fluxtubevolume_',num2str(grid_kind),'_',num2str(coord_system),'_',num2str(NR),'_', ...
       num2str(NPHI),'_',num2str(NTHETA),'.dat'];

    % Set input files
    tetra_grid.TETRA_GRID_NML.n1 = NR;
    tetra_grid.TETRA_GRID_NML.n2 = NPHI;
    tetra_grid.TETRA_GRID_NML.n3 = NTHETA;

    gorilla_applets.GORILLA_APPLETS_NML.filename_fluxtv_load = ['../../../',folder_name_fluxtv,'/',file_name_fluxtv];
    alpha_lifetime.alpha_lifetimenml.filename_alpha_lifetime = ['../results/alpha_life_time_gorilla_poly2_',num2str(NR),'_', ...
                                                                num2str(NPHI),'_',num2str(NTHETA),'.dat'];

    % Write modified Inputfiles
    gorilla.write([gorilla_path,name_run_folder,'/',folder_name,'/gorilla.inp']);
    gorilla_applets.write([gorilla_path,name_run_folder,'/',folder_name,'/gorilla_applets.inp']);
    tetra_grid.write([gorilla_path,name_run_folder,'/',folder_name,'/tetra_grid.inp'])
    alpha_lifetime.write([gorilla_path,name_run_folder,'/',folder_name,'/alpha_lifetime.inp'])

end

  

