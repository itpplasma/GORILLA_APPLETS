% Precomputation of Fluxtube Volumes
% This script precomputes flux tube volumes with an option to only precompute such volumes that are not precomputed yet.

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
  
% read file
gorilla.read();
gorilla_applets.read();
tetra_grid.read();

% General settings for Flux Tube Volume Computation

% Name of path with Fluxtubevolumes
folder_name = 'RUN/fluxtubevolumes';
mkdir(gorilla_path,folder_name);

% Store all names of compued files in cell array
computed_files = struct2cell(dir([gorilla_path,folder_name]));
computed_files = computed_files(1,3:end);


%% General settings for computation of flux tube volume

% Set computation mode to FluxTV computation
gorilla_applets.GORILLA_APPLETS_NML.i_option = 1;

% Choose electros as species
gorilla.GORILLANML.ispecies = 1;

% Grid kind (Stellarator)
grid_kind = 3;
tetra_grid.TETRA_GRID_NML.grid_kind = grid_kind;

% Coordinate system (Symmetry Flux Coordinates)
coord_system = 2;
gorilla.GORILLANML.coord_system = coord_system;


%% Computation of flux tube volumes

% Create RUN folder for fluxtube volume computation
name_run_folder = 'RUN/20220816_fluxtubevolume_computation_02';
mkdir(gorilla_path,name_run_folder);


% Create softlinks
status = system(['ln -s ',gorilla_path,'gorilla_applets_main.x ',gorilla_path,name_run_folder,'/.'],'-echo')
status = system(['ln -s ',gorilla_path,'DATA ',gorilla_path,name_run_folder,'/.'],'-echo')
status = system(['ln -s ',gorilla_path,'DATA/VMEC/wout_23_1900_fix_bdry.nc ',gorilla_path,name_run_folder,'/netcdf_file_for_test.nc'],'-echo')

%grid_vec = [1500];

%grid_vec = [50:50:500];
%grid_vec = [28,32,37,43,49,57,65,75,86,99,114,131,151,174,200];
grid_vec = [28,32];

[NR] = deal(100);

for i = 1:numel(grid_vec)
    
   [NPHI,NTHETA] = deal(grid_vec(i));
   
    % Define filename and path 
    file_name = ['fluxtubevolume_',num2str(grid_kind),'_',num2str(coord_system),'_',num2str(NR),'_', ...
	      num2str(NPHI),'_',num2str(NTHETA),'.dat'];

    % Look, whether FluxTubeVolume was already computed
    boole_precomputed = any(cellfun(@(x) strcmp(x,file_name),computed_files,'UniformOutput',1));
    
    if boole_precomputed
        disp(['The FluxTV configuration ',file_name,' was already computed.'])
    else
        gorilla_applets.GORILLA_APPLETS_NML.filename_fluxtv_precomp = ['../../',folder_name,'/',file_name];
        tetra_grid.TETRA_GRID_NML.n1 = NR;
        tetra_grid.TETRA_GRID_NML.n2 =  NPHI;
        tetra_grid.TETRA_GRID_NML.n3 =  NTHETA;

        % Write modified Inputfiles
        gorilla.write([gorilla_path,name_run_folder,'/gorilla.inp']);
        gorilla_applets.write([gorilla_path,name_run_folder,'/gorilla_applets.inp']);
        tetra_grid.write([gorilla_path,name_run_folder,'/tetra_grid.inp'])
    
        % Run program
        cd([gorilla_path,name_run_folder]);
        status = system('./gorilla_applets_main.x','-echo')
    end
    
end

cd(matlab_path)

  

