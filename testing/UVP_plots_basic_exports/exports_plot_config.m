function options = exports_plot_config
%FUNCTION exports_plot_config
%
%  Syntax:
%    options = exports_plot_config
%
%  Description:
%    Define necessary information to plot this specific project
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
%% Required fields
% Define options structure
options = struct();
options.savefig   = 0;
% Grid via Latitude, longitude, or time
options.grid_type = 'time'; % 'time' 'lat' 'lon'
options.plot_type = 'par'; % 'par' 'zoo' 'ctd'
options.plot_station_data = 0; % 1 = plots individual profiles, 0 = do not plot individual profiles
% Cruise nickname - used for file naming
options.project_name = 'exports';
% Project path used for saving figures and files
fprintf('You will have to change the options.project path to match your system!\n')
keyboard
options.project   = 'D:\MATLAB\UVP_plots_basic_exports\'; 

% Define filename of ZOO and PAR files
%par_file = fullfile(UVP_tools,'testing','p16n_2015_uvpsn009_dataset','export_detailed_20210325_18_53_PAR_odv.txt');
%options.zoo_file = fullfile(UVP_tools,'testing','p16n_2015_uvpsn009_dataset','export_detailed_20210325_18_53_ZOO_odv.txt');
options.par_file = fullfile(options.project,'export_detailed_20210820_16_57_PAR_odv.txt');
options.zoo_file = fullfile(options.project,'export_detailed_20210820_16_57_ZOO_odv.txt');
% generate filename to save matlab formatted UVP data
options.zoo_matfile = fullfile(options.project,[options.project_name '_zoo.mat']);
options.par_matfile = fullfile(options.project,[options.project_name '_par.mat']);

%% Optional fields
% Optional: Cell array of profiles to ignore (must be exact match to "profile"
% read in from ecotaxa file. 
options.exclude_profiles = {'ctd000'};

% Optional: List of fully validated stations for ZOO data plotting
% This will limit zoo data to only the profiles listed in the file
% see formatting requirements in limit_zoo_to_fully_validated_profiles.m
options.validated = fullfile(options.project,'exports_fully_validated_stations.txt');

% Optional: List of sections defined by profile ranges
%** File must be formatted as follows ** 
%"section_p1","section_p2" % first column name, last column name
% section1_p1,section1_p2  % starting profileID, ending profileID
% section2_p1,section2_p2  % starting profileID, ending profileID
% section3_p1,section3_p2  % starting profileID, ending profileID
% section4_p1,section4_p2  % starting profileID, ending profileID
options.filename_sectionsbyprofiles = fullfile(options.project,'sections_by_profile.csv');

% Optional: Bathymetry file
% Must be a structure with fields "lat" "lon" and "z" where "z" is the
% bathymetric depth in meters
%options.bathymetry_file = fullfile(options.project,[options.project_name '_bathymetry.mat']);

% Optional: hardcode which variables you want to plot
% Option I  : Select which fields you want to plot below (see Section 7 of UVP_workflow)
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
options.plots = struct();
options.plots.meansize.title = {'Particle mean size [mm]'};
options.plots.meansize.clims = [0.092167 0.25482];
%% Build files
% Change into project directory so all data will be saved there
if ~exist(options.project,'dir'); mkdir(options.project); end
cd(options.project)

