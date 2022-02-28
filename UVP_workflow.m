function UVP_workflow(toolbox_dir,config_script)
%FUNCTION UVP_WORKFLOW
%
%  Syntax:
%    UVP_workflow(toolbox_dir,config_script)
%
%  Inputs:  
%    toolbox_dir   = directory path for UVP_plots_basic toolbox
%    config_script = configuration script that defines files, and script
%    options
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
%    Code is based on work of Andrew McDonnell, Stephanie O'Daly,
%    and Jessica Pretty. 
% 
%%
% THIS IS A TEST TO SEE IF CHANGES SHOW UP ON AMCDONNELL-LAB ACCOUNT
close all
fprintf(' ********************************************************\n')
fprintf(' ********************************************************\n')
fprintf(' ADD PLOTS:\n')
fprintf('     1. diurnal patterns day/night +/- 1 hour sunset\n')
fprintf('     2. Pie charts of zoo living/non-living\n')
fprintf('     3. euphotic vs non-euphotic\n')
fprintf('     4. break into futher depth sections?\n')
fprintf('     5. Others from EXPORTS UVP Results webinar Aug3-2020\n') %https://docs.google.com/presentation/d/1tHKvPtbA6SyDpSdu59G_l5s-BkLOG_uI/edit?usp=sharing&ouid=110344274478839061352&rtpof=true&sd=true
fprintf(' ********************************************************\n')

%% 0 | Make sure you're in the UVP_tools repository
% *** CHANGE THIS DIRECTORY PATH TO MATCH YOUR LOCAL COMPUTER ***
addpath(genpath(toolbox_dir));  % Path to UVP_plots_basic/utility & necessary other folders that contains dependent functions

%% 1 | Read configuration script
% first, parse the path to the configure script
[project_folder,script_name,~] = fileparts(config_script);

try
  cd(project_folder) % changed directory into project folder
  eval(['options = ' script_name])% evaulate the script to read in info
catch
  fprintf('Failed trying to run configure script: %s \n',config_script)
  fprintf('Check the path is correct\n');
end

% Check that necessary components are there
if ~isfield(options,'savefig')
  options.savefig   = 0; % default to NO
end
if ~isfield(options,'grid_type')
  options.grid_type = 'time'; % 'time' 'lat' 'lon'
end
if ~isfield(options,'plot_type')
  options.plot_type = 'par'; % 'par' 'zoo' 'ctd'
end
if ~isfield(options,'plot_station_data')
  options.plot_station_data = 0; % default to NO
end

%% 2 | Read PAR file and calculate parameters from particle abundance and biovolume data
if exist(options.par_matfile,'file')
  load(options.par_matfile);
else
  % 1a| Read ODV PAR file exported from Ecotaxa
  [par_table, par_info] = UVP_read_odv_ecotaxa_exported_par(options.par_file);
  % 1b| Calculate more parameters from particle abundance and biovolume data
  [par,par_info] = UVP_calculate_PAR_fields(par_table,par_info);

  % 1c| Save data to file, the -v7.3 is just a compression flag
  save(options.par_matfile,'par','par_info','-v7.3');
end

%% 3 | Read ZOO file with option to limit plotting to fully validated stations
if exist(options.zoo_matfile,'file')
  fprintf('Loading ZOO data %s\n',options.zoo_matfile);
  load(options.zoo_matfile);
  
else
  % 2a| Read ZOO file exported from Ecotaxa
  [zoo, zoo_info] = UVP_read_odv_ecotaxa_exported_zoo(options.zoo_file);
  
  % 2b| Save zoo data to file, the -v7.3 is just a compression flag
  save(options.zoo_matfile,'zoo','zoo_info','-v7.3');
end
% Limit to fully validated stations
% Read in file that contains a list of fully validated profile names
if exist(options.validated,'file')
  zoo = limit_zoo_to_fully_validated_profiles(options.validated,zoo);
end
%  add project for filename creation and clarify
try zoo.header = options.project_name;end
try par.header = options.project_name;end

%% 4 | Read CTD data
% Still need to incorporate this...
% Read_hydrographic_cruise_hyfiles

%% 5 | Remove unwanted profiles
if isfield(options,'exclude_profiles')
  remove_profiles_par = ismember(par.profile,options.exclude_profiles);
  remove_profiles_zoo = ismember(zoo.profile,options.exclude_profiles);
  % Convert to table to remove all profile data
  if isfield(par,'header'); par = rmfield(par,'header'); end
  if isfield(zoo,'header'); zoo = rmfield(zoo,'header'); end
  par = struct2table(par);
  zoo = struct2table(zoo);
  % remove unwanted profile data
  par(remove_profiles_par,:) = [];
  zoo(remove_profiles_zoo,:) = [];
  % convert back to structure
  par = table2struct(par,'ToScalar',true);
  zoo = table2struct(zoo,'ToScalar',true);
end

%% 7 | Set data to ZOO OR PAR OR CTD data
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
keyboard
%% 7 | Manually select OR Hardcode which variables you want to plot
% Option I  : Select which fields you want to plot below (see Section 7)
%   OR
% Option II : hardcode using the following format in the config_script
%   "plots" is a structure where fieldnames are the exact fieldnames of the
%   variables you want to plot. In addition, each field has "title" and
%   "clims" fields.
% EXAMPLE #1: data.tot_par_abundance (i.e. "total particle abundance")
%   options.tot_par_abundance       = struct();
%   options.tot_par_abundance.title = {'Abundance (102-203um)[#/L]'};
%   options.tot_par_abundance.clims = [5.35 57.86];
% EXAMPLE 2: data.VSD (i.e. "particle biovolume"). Since VSD is an NxM
%   array, you need to define which sizes you want to use.
%   options.VSD       = struct();
%   options.VSD.title = {'Biovolume (102-128_um)[ppm]'; 'Biovolume (128-161_um)[ppm]'; 'Biovolume (161-203_um)[ppm]'};
%   options.VSD.clims = [[0.0002 0.00366];[0.000559 0.0094];[0.00058 0.0106]];
% EXAMPLE 3: ZOO data plotting abundance and biovolume for select taxa
%   plots = struct();
%   options.abund.title = {'abund_Rhizaria_Harosa';'abund_Crustacea_Arthropoda'};
%   options.abund.clims = [[0 20];[0 20]];
%   options.biovol.title = {'biovol_Rhizaria_Harosa';'biovol_Crustacea_Arthropoda'};
%   options.biovol.clims = [[0 20];[0 20]];

%% 8 | Determine which variables to plot
if ~isfield(options,'plots') || isempty(options.plots)
  try
    options.plots = UVP_select_paramets_to_plot(data,data_info);
  catch
    fprintf('Could not run UVP_select_paramets_to_plot.m script\n')
    fprintf('so... you will need to define the variables manually, or debug\n')
  end
end

%% 9 | Remove two deepest depth bins in every profile
% default to removing last two depth bins
[data] = uvp_remove_last_two_depth_bins(data);

% Add basic data fields to options structure
options.profile = data.profile;
options.time    = data.datenum;
options.lat     = data.latitude;
options.lon     = data.longitude;
options.depth   = data.Depth;

%% 10 | Decide if want to plot sections of data, or the whole range
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

%% 11 | Load bathymetry
try
  if isfield(options,'bathymetry_file') && exist(options.bathymetry_file,'file')
    load(options.bathymetry_file,'bathy')
    % Interpolate bathymetry to cruise latitude longitude
    data.bathy = interp2(bathy.lon,bathy.lat,bathy.z,data.longitude,data.latitude,'linear','extrap');
  else
    % bathy = structure with latitude, longitude, and bathymetry positions
    bathy = get_bathymetry(data.longitude,data.latitude,options.project_name);
    data.bathy = interp2(bathy.lon,bathy.lat,bathy.z,data.longitude,data.latitude,'linear','extrap');
  end
  % Update options
  if isfield(data,'bathy')
    options.bathy = data.bathy;
  end
catch
  fprintf('Could not load or download bathymetry\n')
  fprintf('...see instructions on github page for details if you want to include bathymetry\n')
end

%% 12 | PLOTS | Standard PAR 2D vertical profile plots
%% 12a | Plot standard 3 panel total particle abundance, total particle biovolume, and slope of PSD
if strcmp(options.plot_type,'par')
  plot_uvp_multipanel(par,par_info,{'tot_par_abundance' 'tot_par_biovolume' 'slope_b'},options);
end

%% 12b | Plot vertical profiles of different sizes at each station
% NSD is short for number size distribution and is in units of #/L
if options.plot_station_data
  plot_uvp_NSD(par,par_info,options);
end

%% 13 | PLOTS | Generate plots for each field
% This just loops through selected fields
% Comment out specific plotting functions if you want to skip them
% For input requirements for functions, see the function documentation
fields_to_plot = fieldnames(options.plots);
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
    for nn = 1:size(options.plots.(field).title,1)
      options.plot_title = options.plots.(field).title{nn};   % Pull out title
      options.data_clims = options.plots.(field).clims(nn,:); % Pull out color limits
      options.data       = data.(field)(:,nn);
      fprintf('\n-----------------------------------------\n')
      fprintf('Plotting %s: %s\n',field,options.plot_title)
      fprintf('-----------------------------------------\n')
      %% PLOTS OF NON-GRIDDED DATA
      if options.plot_station_data
        % 1 | Waterfall plot - shows evolution of profiles throughout cruise
        plot_data_vs_depth_waterfall(options);
        
        % 2 | Plot of variable vs depth
        plot_data_vs_depth(options);
      end
      
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
      LAT = griddata(x,options.depth(xrng),options.lat(xrng),X,Z); % griddata okay here because nothing fancy, just latitude
      LON = griddata(x,options.depth(xrng),options.lon(xrng),X,Z); % griddata okay here because nothing fancy, just longitude
      
      % 2 | 2D confour plot
      fig = plot_2d_gridded_transect(X,Z,DAT,options);
      plot_map_inset(fig,options);
      
      % 3 | 3D contour plot
      %fig = plot_3d_gridded_transect(LON,LAT,Z,DAT,options,bathy);
      %plot_map_inset(fig,options);
      
    end  %% LOOP THROUGH NxM DATA
    % ---------------------------------------------------------------------
  end %% LOOP THROUGH SECTIONS OF DATA
  % -----------------------------------------------------------------------
end %% LOOP THROUGH FIELDS TO PLOT



