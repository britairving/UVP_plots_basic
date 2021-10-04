function UVP_workflow
%FUNCTION UVP_WORKFLOW
%
%  Syntax:
%    UVP_WORKFLOW
%
%  Description:
%    UVP_workflow is the wrapper script that you can set up to call
%    plotting functions.
%
%  Citation:
%    Brita Irving, (2021), UVP_plots_basic, GitHub repository,
%    https://github.com/britairving/UVP_plots_basic
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%
close all

%% 0 | Make sure you're in the UVP_tools repository
% *** CHANGE THIS DIRECTORY PATH TO MATCH YOUR LOCAL COMPUTER ***
if ismac
  UVP_tools = '/Users/bkirving/Documents/MATLAB/UVP_plots_basic'; 
elseif ispc
  UVP_tools = 'D:\MATLAB\UVP_plots_basic';
end
addpath(genpath(UVP_tools));  % Path to UVP_tools/utility & necessary other folders that contains dependent functions
cd(UVP_tools)
%% 1 | Define project nickname & Ecotaxa filenames
% UVP data must be exported in detailed ODV format
% < https://ecotaxa.obs-vlfr.fr/part/ >
% Define options structure
options.savefig   = 0;
options.project   = fullfile(pwd,'testing','exports'); % Project nickname use for saving figures and files
% Grid via Latitude, longitude, or time
options.grid_type = 'lat'; % 'time' 'lat' 'lon'
options.plot_type = 'par'; % 'par' 'zoo' 'ctd'
% List of fully validated stations for ZOO data plotting
options.validated = [options.project '_fully_validated_stations.txt'];
options.filename_sectionsbyprofiles = fullfile(options.project,'sections_by_profile.csv');
% Define filename of ZOO and PAR files
%par_file = fullfile(UVP_tools,'testing','p16n_2015_uvpsn009_dataset','export_detailed_20210325_18_53_PAR_odv.txt');
%zoo_file = fullfile(UVP_tools,'testing','p16n_2015_uvpsn009_dataset','export_detailed_20210325_18_53_ZOO_odv.txt');
par_file = fullfile(UVP_tools,'testing','exports','export_detailed_20210820_16_57_PAR_odv.txt');
zoo_file = fullfile(UVP_tools,'testing','exports','export_detailed_20210820_16_57_ZOO_odv.txt');

% Change into project directory so all data will be saved there
if ~exist(options.project,'dir'); mkdir(options.project); end
cd(options.project)

% generate filename to save matlab formatted UVP data
zoo_matfile = [options.project '_zoo.mat'];
par_matfile = [options.project '_par.mat'];

%% 2 | Read PAR file and calculate parameters from particle abundance and biovolume data
if exist(par_matfile,'file')
  load(par_matfile);
else
  % 1a| Read ODV PAR file exported from Ecotaxa
  [par_table, par_info] = UVP_read_odv_ecotaxa_exported_par(par_file);
  % 1b| Calculate more parameters from particle abundance and biovolume data
  [par,par_info] = UVP_calculate_PAR_fields(par_table,par_info);
  % 1c| add project for filename creation and clarify
  par.header = options.project;
  % 1d| Save data to file, the -v7.3 is just a compression flag
  save(par_matfile,'par','par_info','-v7.3');
end

%% 3 | Read ZOO file with option to limit plotting to fully validated stations
if exist(zoo_matfile,'file')
  fprintf('Loading ZOO data %s\n',zoo_matfile);
  load(zoo_matfile);
  
else
  % 2a| Read ZOO file exported from Ecotaxa
  [zoo, zoo_info] = UVP_read_odv_ecotaxa_exported_zoo(zoo_file);
  
  % 2b| add project for filename creation and clarify
  zoo.header = options.project;
  
  % 2c| Save zoo data to file, the -v7.3 is just a compression flag
  save(zoo_matfile,'zoo','zoo_info','-v7.3');
end
% Limit to fully validated stations
% Read in file that contains a list of fully validated profile names
if exist(options.validated,'file')
  zoo = limit_zoo_to_fully_validated_profiles(options.validated,zoo);
end

%% 4 | Read CTD data
% Still need to incorporate this...
% Read_hydrographic_cruise_hyfiles
 
%% 5 | Set data to ZOO OR PAR OR CTD data
% If you want to plot ZOO variables, data = zoo; etc
switch options.plot_type
  case 'zoo'
    data      = zoo;      % or zoo, or ctd
    data_info = zoo_info; % or zoo_info, or ctd_info'
    data.abund(isnan(data.abund))  = 0;
    data.biovol(isnan(data.abund)) = 0;
  case 'par'
    data      = par;      % or zoo, or ctd
    data_info = par_info; % or zoo_info, or ctd_info
  case 'ctd'
    fprintf('CTD not set up yet..\n')
    keyboard
    data      = ctd;      % or zoo, or ctd
    data_info = ctd_info; % or zoo_info, or ctd_info
  otherwise
    error('plot_type not recognized: choices are zoo, par, or ctd')
end

%% 6 | Manually select OR Hardcode which variables you want to plot
% Option I  : Select which fields you want to plot below (see Section 7)
%   OR
% Option II : hardcode using the following format
%   "plots" is a structure where fieldnames are the exact fieldnames of the
%   variables you want to plot. In addition, each field has "title" and
%   "clims" fields.
% EXAMPLE #1: data.tot_par_abundance (i.e. "total particle abundance")
%   plots.tot_par_abundance       = struct();
%   plots.tot_par_abundance.title = {'Abundance (102-203um)[#/L]'};
%   plots.tot_par_abundance.clims = [5.35 57.86];
% EXAMPLE 2: data.VSD (i.e. "particle biovolume"). Since VSD is an NxM
%   array, you need to define which sizes you want to use.
%   plots.VSD       = struct();
%   plots.VSD.title = {'Biovolume (102-128_um)[ppm]'; 'Biovolume (128-161_um)[ppm]'; 'Biovolume (161-203_um)[ppm]'};
%   plots.VSD.clims = [[0.0002 0.00366];[0.000559 0.0094];[0.00058 0.0106]];
% EXAMPLE 3: ZOO data plotting abundance and biovolume for select taxa
%   plots = struct();
%   plots.abund.title = {'abund_Rhizaria_Harosa';'abund_Crustacea_Arthropoda'};
%   plots.abund.clims = [[0 20];[0 20]];
%   plots.biovol.title = {'biovol_Rhizaria_Harosa';'biovol_Crustacea_Arthropoda'};
%   plots.biovol.clims = [[0 20];[0 20]];

plots = struct();
plots.meansize.title = {'Particle mean size [mm]'};
plots.meansize.clims = [[0.092167 0.25482]];
%% 7 | Determine which variables to plot
if ~exist('plots','var') || isempty(plots)
  try
    plots = UVP_select_paramets_to_plot(data,data_info);
  catch
    fprintf('Could not run UVP_select_paramets_to_plot.m script\n')
    fprintf('so... you will need to define the variables manually, or debug\n')
  end
end


%% 8 | Remove two deepest depth bins in every profile
% default to removing last two depth bins
[data] = uvp_remove_last_two_depth_bins(data);

% Add basic data fields to options structure
options.profile = data.profile;
options.time    = data.datenum;
options.lat     = data.latitude;
options.lon     = data.longitude;
options.depth   = data.Depth;


%% 9 | Decide if want to plot sections of data, or the whole range
try
  fprintf('\nChoose how to visualize data\n')
  fprintf(' <1> All data (Default)\n')
  fprintf(' <2> Split into sections\n')
  section_chc = input('Enter choice: ');
  if isempty(section_chc) || section_chc == 1
    num_sections = 1;
  else % CHOOSE
    options.sections = select_data_sections(options);
    % Returns NxM matrix of [section#,indices]
    num_sections = size(options.sections,1);
  end
catch
  num_sections = 1;
end

%% 10 | Load bathymetry
try
  if exist([data.header '_bathymetry.mat'],'file')
    load([data.header '_bathymetry.mat'])
    if isfield(bathy,'bathy_transect')
      data.bathy = interp1(bathy.lat_transect,bathy.bathy_transect,data.latitude);
    end % check if available
  else
    % bathy = structure with latitude, longitude, and bathymetry positions
    bathy = get_bathymetry(data.longitude,data.latitude,data.header);
    data.bathy = interp1(bathy.lat_transect,bathy.bathy_transect,data.latitude);
  end
  % Update options
  if isfield(data,'bathy')
    options.bathy = data.bathy;
  end
catch
  fprintf('Could not load or download bathymetry\n')
  fprintf('...see instructions on github page for details if you want to include bathymetry\n')
end

%% 11 | PLOTS | Standard PAR 2D vertical profile plots
%% 11a | Plot standard 3 panel total particle abundance, total particle biovolume, and slope of PSD
%plot_uvp_multipanel(par,par_info,{'tot_par_abundance' 'tot_par_biovolume' 'slope_b'},options);


%% 11b | Plot vertical profiles of different sizes at each station
% NSD is short for number size distribution and is in units of #/L
%plot_uvp_NSD(par,par_info,options);

%% 12 | PLOTS | Generate plots for each field
% This just loops through selected fields
% Comment out specific plotting functions if you want to skip them
% For input requirements for functions, see the function documentation
fields_to_plot = fieldnames(plots);
for nfield = 1:numel(fields_to_plot)
  field = fields_to_plot{nfield};
  if ~isfield(data,field)
    fprintf('%s is not available in data\n',field)
    continue
  end
  % -----------------------------------------------------------------------
  %% Loop through sections of data
  for nrng = 1:num_sections
    % Pull out indices of current section
    if num_sections > 1
      xrng = options.sections(nrng,1):options.sections(nrng,2);
      options.section_idx = xrng; % used in plot_map_inset.m
    else
      xrng = 1:numel(options.time);
    end
    % ---------------------------------------------------------------------
    %% Loop through available data in case it is NxM array
    % This is the case for VSD, NSD, DNSD, DVSD, & ZOO parameters.
    for nn = 1:size(plots.(field).title,1)
      options.plot_title = plots.(field).title{nn};   % Pull out title
      options.data_clims = plots.(field).clims(nn,:); % Pull out color limits
      options.data       = data.(field)(:,nn);
      fprintf('\n-----------------------------------------\n')
      fprintf('Plotting %s: %s\n',field,options.plot_title)
      fprintf('-----------------------------------------\n')
      %% PLOTS OF NON-GRIDDED DATA
      % 1 | Waterfall plot - shows evolution of profiles throughout cruise
      %plot_data_vs_depth_waterfall(options);

      % 2 | Plot of variable vs depth
      %plot_data_vs_depth(options);
      
      
      %% --------------------------------------------------------------------
      %% PLOTS OF GRIDDED DATA
      % 1 | Grid data as NxM - X vs Y
      % where Y = Depth and X = latitude, longitude, or time
      switch options.grid_type
        case 'time'; x = data.datenum(xrng);
        case 'lat';  x = data.latitude(xrng);
        case 'lon';  x = data.longitude(xrng);
        otherwise
          fprintf('unexpected grid type: %s\n',grid_type)
          fprintf('grid type should be "time" "lat" or "lon"\n')
          keyboard
      end
      
      [X,Z,DAT] = UVP_griddata(x,options.depth(xrng),options.data(xrng));
      LAT = griddata(x,options.depth(xrng),options.lat(xrng) ,X,Z); % griddata okay here because nothing fancy, just latitude
      LON = griddata(x,options.depth(xrng),options.lon(xrng),X,Z); % griddata okay here because nothing fancy, just longitude
      
      % 2 | Waterfall plot - shows evolution of data throughout cruise
      plot_gridded_data_vs_depth_waterfall(X,Z,DAT,options);
      keyboard
      % 3 | 2D confour plot
      fig = plot_2d_gridded_transect(X,Z,DAT,options);
      plot_map_inset(fig,options);
      keyboard
      % 4 | 3D contour plot
      fig = plot_3d_gridded_transect(LON,LAT,Z,DAT,options,bathy);
      plot_map_inset(fig,options);
      keyboard
    end  %% LOOP THROUGH NxM DATA
    % ---------------------------------------------------------------------
  end %% LOOP THROUGH SECTIONS OF DATA
  % -----------------------------------------------------------------------
end %% LOOP THROUGH FIELDS TO PLOT



