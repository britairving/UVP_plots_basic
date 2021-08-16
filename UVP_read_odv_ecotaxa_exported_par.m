function [odv, info] = UVP_read_odv_ecotaxa_exported_par(par_file)
%% function UVP_read_odv_ecotaxa_exported_par
%
%  Syntax:
%    summary = UVP_read_odv_ecotaxa_exported_par(par_file)
%
%  Description:
%    Reads detailed ODV file exported from Ecotaxa https://ecotaxa.obs-vlfr.fr/
%
%  Inputs:
%    par_file = Path to detailed ODV file exported from Ecotaxa
%
%  Outputs:
%    odv = table containing data exported from Ecotaxa
%
%  Example:
%
%  Dependencies:
%   Matlab R2016b or later (uses detectImportOptions.m and contains.m)
%
%  Author: Brita Irving <bkirving@alaska.edu>
%% 0 | Set script flags
% See section 6 | DO NOT DELETE higher size bins after first size bin with real data
% delete_all_empty_size_bins == 1 | Deletes ALL empty size bins
% delete_all_empty_size_bins == 0 | Keeps all size bins after the first
% one that contains data, even if that size bin contains no real data. E.g.
% could be useful for provenance or data archival purposes. All size bins
% below the first one to contain data are deleted.
delete_all_empty_size_bins = 1; 


%% 0 | See if par_file exists
% File must have been exported from Ecotaxa Particle Module in Detailed ODV
% format.  <https://ecotaxa.obs-vlfr.fr/part/>
% IF file does not exist, or was not passed through as an input argument
if nargin < 1
  error('Must pass through full file path to detailed ODV file exported from Ecotaxa')
elseif ~exist(par_file,'file')
  fprintf(' ** File not found: %s **\n',par_file)
  fprintf('File must have been exported from Ecotaxa Particle Module in Detailed ODV format\n')
  error('Could not find par file. Expected full file path to detailed ODV file exported from Ecotaxa')
end

%% 1 | Read column names and header information from ODV file
fprintf('\n  Reading ODV Particle data from... %s\n',par_file)
% First - get table configuration using Matlab builtin detectImportOptions
% E.g. delimiter can be ';' or '\t' (tab)
table_options = detectImportOptions(par_file);
% Read header
% column names usually at line 7, but go to line 20 in case extra
% header information is in file
hdr = cell(20,1);
fileID = fopen(par_file);
for nline = 1:20
  hdr{nline} = fgetl(fileID);
end
fclose(fileID);
column_hdrline = find(contains(hdr,'Cruise'));
% resize header
hdr = hdr(1:column_hdrline);
% the "cols_original" cell contains the exact column names 
cols_original = strsplit(hdr{column_hdrline},table_options.Delimiter);  % split hdr by odv delimiter
cols = cols_original; 

for ic = 1:numel(cols)
  cols{ic} = strrep(cols{ic},' ','_');
end

% ecotaxa changed format of column names so need to adapt here
if any( contains(cols,'_[#_l-1]') )
  cols = strrep(cols,'_[#_l-1]','[#/L]');
end
if any (contains(cols,')_[mm3_l-1]') )
  cols = strrep(cols,')_[mm3_l-1]',')[ppm]'); % [mm^3/L] = [ppm]
end
% Replace micro character with "u", micro character is char(181) in matlab
cols = strrep(cols,char(181),'u');
%% 2 | Reformat variable names
meta_idx = find(contains(table_options.VariableNames,'METAVAR'));
% Remove unnecessary column name text
for nm = 1:numel(meta_idx)
  irm_string = strfind(table_options.VariableNames{meta_idx(nm)},'_METAVAR');
  table_options.VariableNames{meta_idx(nm)} = table_options.VariableNames{meta_idx(nm)}(1:irm_string-1);
  cols{meta_idx(nm)} = cols{meta_idx(nm)}(1:irm_string-1);
  if contains(table_options.VariableNames{meta_idx(nm)},'DOUBLE')
    table_options.VariableTypes{meta_idx(nm)} = 'double';
  end
end
% make sure all variables read in correctly as number and not strings
table_options.VariableTypes(meta_idx(end):end) = {'double'};
% remove other primary variable identifier - usually on Depth [m] variable
table_options.VariableNames = erase(table_options.VariableNames,'__PRIMARYVAR_DOUBLE');
cols = erase(cols,':PRIMARYVAR:DOUBLE');
table_options.VariableNames(1:meta_idx(end)) = lower(table_options.VariableNames(1:meta_idx(end)));
table_options.VariableNames{contains(table_options.VariableNames,'yyyy_mm_ddhh_mm')} = 'datetime';

% Changing these variable names is not 100% necessary, but  keeping it for
% now just for provenance and backward compatability. 
table_options.VariableNames = erase(table_options.VariableNames,{'_degrees_north_' '_degrees_east_'}); % Remove units on latitude and longitude
table_options.VariableNames{contains(table_options.VariableNames,'station')} = 'profile';
table_options.VariableNames{contains(table_options.VariableNames,'Depth_m')} = 'Depth';
table_options.VariableNames{contains(table_options.VariableNames,'SampledVolume_L_')} = 'SampledVolume';

%% 3 | Read ODV particle data
try
  odv = readtable(par_file,table_options);
catch
  fprintf('could not read %s with table_options\n',par_file);
  fprintf('stopping here in keyboard mode so can troubleshoot..\n')
  keyboard
end

%% 4 | Delete all "extrames" variables
rm_extrames = contains(odv.Properties.VariableNames,'extrames');
odv(:,rm_extrames) = [];
cols(rm_extrames)  = [];

%% 5 | Fill meta variables
% ODV format has a single row of metadata per station. 
fnames = fieldnames(odv);
for iv = 1:numel(meta_idx)
  vname = fnames{meta_idx(iv)};
  odv.(vname) = fillmissing(odv.(vname),'previous');
end % loop through meta variables

%% 6 | Identify all empty variables - i.e. no real values
% This can include size bins, ctd variables, extrames variables, etc. Any
% column that contains data no finite data is flagged for removal.
fnames = fieldnames(odv);
rm_fields = [];
for irm = meta_idx(end)+1:size(odv,2)
  nam = fnames{irm};
  if ~iscell(odv.(nam))
    if all(isnan(odv.(nam)))
      rm_fields = [rm_fields; irm];
    elseif all(odv.(nam) == 0)
      rm_fields = [rm_fields; irm];
    end
  end
end

%% 7 | DO NOT DELETE higher size bins after first size bin with real data
if delete_all_empty_size_bins == 0
  % For example, in many UVP5 HD datasets, the 50.8-64um contains data, but
  % the next size bin does not. For the non-HD units, the 102-128um size bin
  % contains data but the 128-161um does not.
  % Detection limit for UVP is approximately 100 um but there is commonly
  % data in the 50.8-64um size bin and we want to maintain all columns with
  % data for posterity sake.
  nperL  = find(contains(cols,'[#/L]')); % Indices of number concentration variables
  biovol = find(contains(cols,'[ppm]')); % Indices of biovolume variables
  % Find first size class that contains data
  nperL_data     = table2array(odv(:,nperL));
  idx_firstgood  = find(any(isfinite(nperL_data),1),1);
  idx_keep_nperL = nperL(idx_firstgood):nperL(end);
  % Find first size class that contains data
  biovol_data    = table2array(odv(:,biovol));
  % Biovolume fields with no value are exported from ecotaxa as 0, instead of
  % blank. Therefore, replace all biovolume == 0 with NaN to avoid confusion.
  biovol_data(biovol_data == 0) = NaN;
  idx_firstgood   = find(any(isfinite(biovol_data),1),1);
  idx_keep_biovol = biovol(idx_firstgood):biovol(end);
  ignore_empty    = ismember(rm_fields, [idx_keep_nperL idx_keep_biovol] );
  rm_fields(ignore_empty) = []; % keep the columns meeting the above criteria
end

%% 7 | Remove empty fields 
odv(:,rm_fields) = [];
cols(rm_fields) = [];

%% 8 | Remove fields named "extrames" because those could be anything
rm_fields = find(contains(cols,'extrames'));
odv(:,rm_fields) = [];
cols(rm_fields) = [];

%% 9 | Parse size bin limits from fieldnames
% split into categories #/L & biovolume
nperL  = find(contains(cols,'[#/L]'));
biovol = find(contains(cols,'[ppm]'));

% Pull out size class limits and theoretical bin centers
num_sizes = numel(nperL);
vsize = struct();
vsize.strnam = cell(num_sizes,1);   % string containing only bin limits
vsize.NSDnam = cell(num_sizes,1);   % simple field name for abundance variables
vsize.VSDnam = cell(num_sizes,1);   % simple field name for biovolume variables
vsize.sizemm = nan(num_sizes,2);    % bin range in mm
vsize.sizeum = nan(num_sizes,2);    % bin range in microns
vsize.sizeav_mm = nan(num_sizes,1); % bin center in mm
vsize.sizeav_um = nan(num_sizes,1); % bin center in microns
for in = 1:num_sizes
  % isolate string containing only bin limits
  vsize.strnam{in} = erase(cols{nperL(in)},{'(', ')','LPM_','[#/L]','_'});
  % generate basic field to use later
  vsize.NSDnam{in} = ['LPM' num2str(in)];             % simple field name for abundance variables
  vsize.VSDnam{in} = ['LPM_biovolume' num2str(in)];   % simple field name for biovolume variables
  % pull out sizes from variable name
  if contains(vsize.strnam{in} ,'>') % > 26mm
    str_size1 = erase(vsize.strnam{in},{'mm', '>'});
    sz1 = str2double(str_size1);
    vsize.sizemm(in,:) = [sz1 NaN];
  else
    str_sizes = strsplit(vsize.strnam{in},'-');
    str_size1 = erase(str_sizes{1},{'mm', 'um'});
    str_size2 = erase(str_sizes{2},{'mm', 'um'});
    sz1 = str2double(str_size1);
    sz2 = str2double(str_size2);
    % convert numbers to mm
    if contains(vsize.strnam{in},'um') && contains(vsize.strnam{in},'mm') % only first number in [um]
      vsize.sizemm(in,:) = [sz1./1000 sz2];
    elseif contains(vsize.strnam{in},'um') % both numbers in [um]
      vsize.sizemm(in,:) = [sz1./1000 sz2./1000];
    elseif contains(vsize.strnam{in},'mm') % both numbers in [mm]
      vsize.sizemm(in,:) = [sz1 sz2];
    else
      fprintf('unexpected scenario.. stopping to check here\n')
      keyboard
    end
  end
  % Calculate size range in microns
  vsize.sizeum(in,:) = vsize.sizemm(in,:).*1000; % convert mm to microns
  % Calculate bin center
  % bin_center = 10_to_the_power_of_the_mean_of_the_log10_transformed_size_bin_limits_in_micrometers
  % this is really just a theoretical bin center
  % Pull out log10 center of the bin limits in mm
  meansz_mm = logspace(log10(vsize.sizemm(in,1)),log10(vsize.sizemm(in,2)),3); % returns limit1, center, limit2
  vsize.sizeav_mm(in) = meansz_mm(2); % pull out center of size range
  % Pull out log10 center of the bin limits in microns
  meansz_um = logspace(log10(vsize.sizeum(in,1)),log10(vsize.sizeum(in,2)),3); % returns limit1, center, limit2
  vsize.sizeav_um(in) = meansz_um(2); % pull out center of size range
end

% Pull out size bin to include in human readable description
size_bin_limits = strrep(vsize.strnam,'um', ' micrometers'); 
size_bin_limits = strrep(size_bin_limits,'mm', ' millimeters');
abun_desc   = 'Abundance of particles with an equivalent spherical diameter between';
biovol_desc = 'Biovolume of particles with an equivalent spherical diameter between'; % 'Biovolume, or volume size distribution, for particles';
abun_desc   = strcat(abun_desc,  {' '},size_bin_limits)';
biovol_desc = strcat(biovol_desc,{' '},size_bin_limits)';

%% 10 | Simplify fieldnames
% Currently, the column names in the table are generated automaticaly from
% matlab and are not very easily human readable
% For example: "LPM_50_8_64_m___L_1_" == "LPM (50.8-64_um) [#/L]"
nperL  = find(contains(cols,'[#/L]')); % Indices of number concentration variables
biovol = find(contains(cols,'[ppm]')); % Indices of biovolume variables
for ii = 1:numel(nperL)
  odv.Properties.VariableNames{nperL(ii)}  = vsize.NSDnam{ii};
  odv.Properties.VariableNames{biovol(ii)} = vsize.VSDnam{ii};
end

%% 11 | Add metadata for each variable in a structure
descs  = {'Cruise' 'Station' 'Profile' 'Dataowner' 'UVP_Rawfilename' 'Instrument' 'SerialNumber' 'CTDrosettefilename' 'Datetime_UTC'        'Latitude'      'Longitude'   'Depth' 'Sampled_Volume'}; 
units  = {'none'   'none'    'none'    'none'      'none'            'none'       'none'         'none'               'dd-mm-yyyy HH:MM:SS' 'degrees_north' 'degrees_east' 'm'    'L'};
descs  = [descs abun_desc biovol_desc];
units  = [units repmat({'#/L'},size(nperL)) repmat({'mm^3/L'},size(biovol))];

%% 12 | Add metadata for CTD fields (if available)
% Ecotaxa requires CTD to be submitted in a very specific format
% The list of pre-defined column names is (at time of script creation) :
%     Chloro Fluo [mg Chl/m3]
%     Conductivity [mS/cm]
%     Cpar [%]
%     Depth [salt water, m]
%     Fcdom factory [ppb QSE]
%     In situ Density Anomaly [kg/m3]
%     Neutral Density [kg/m3]
%     Nitrate [umol/l]
%     Oxygen [umol/kg]
%     Oxygen [ml/l]
%     PAR [umol m-2 s-1]
%     Part backscattering coef 470 nm [m-1]
%     Pot. Temperature [degC] (any ref.)
%     Potential Density Anomaly [kg/m3]
%     Potential Temperature [degC]
%     Practical Salinity [psu]
%     Practical Salinity from Conductivity
%     Pressure in Water Column [db] => MANDATORY for DEPTH PROFILES
%     QC Flag
%     Sound Speed c [m/s]
%     SPAR [umol m-2 s-1]
%     Temperature [degC]
%     Time [yyyymmddhhmmssmmm] => MANDATORY FOR TIME DISPLAY

idx_lastbiovolume = find(contains(cols,'biovolume'),1,'last');
if numel(cols) > idx_lastbiovolume
  % Loop through extra fields and build information 
  ctd_indices = idx_lastbiovolume+1:numel(cols);
  for nidx = 1:numel(ctd_indices)
    switch cols{ctd_indices(nidx)}
      case 'chloro_fluo_[mg_chl_m-3]'
        odv.Properties.VariableNames{ctd_indices(nidx)} = 'chlorophyll';
        descs = [descs 'Chlorophyll'];
        units = [units 'mg Chl/m^3'];
      case 'oxygen_[umol_kg-1]'
        odv.Properties.VariableNames{ctd_indices(nidx)} = 'oxygen';
        descs = [descs 'Oxygen'];
        units = [units 'umol/kg'];
      case 'practical_salinity_[psu]'
        odv.Properties.VariableNames{ctd_indices(nidx)} = 'salinity';
        descs = [descs 'Salinity'];
        units = [units 'psu'];
      case 'pressure_[db]'
        odv.Properties.VariableNames{ctd_indices(nidx)} = 'pressure';
        descs = [descs 'Pressure'];
        units = [units 'dbar'];
      case 'temperature_[degc]'    
        odv.Properties.VariableNames{ctd_indices(nidx)} = 'temperature';
        descs = [descs 'Temperature'];
        units = [units 'C'];
      otherwise
        fprintf('Need to add description and units of "%s"\n',cols{ctd_indices(nidx)})
        keyboard
    end
  end % LOOP THROUGH CTD FIELDS
end % IF THERE ARE CTD FIELDS

%% 13 | Generate structure to store detailed variable information
info =  struct();
info.size_bins = vsize;
for nfield = 1:numel(odv.Properties.VariableNames)
  sfield = odv.Properties.VariableNames{nfield};
  info.(sfield).name      = cols{nfield};
  info.(sfield).unit      = units{nfield};
  info.(sfield).long_name = descs{nfield};
  if ismember(sfield,vsize.NSDnam)
    idx_match = strcmp(vsize.NSDnam,sfield);
    info.(sfield).size_bin_st  = vsize.strnam{idx_match}; % string containing size bin limits
    info.(sfield).size_bin_mm  = vsize.sizemm(idx_match,:); % Size bin limits in mm
    info.(sfield).size_mid_mm  = vsize.sizeav_mm(idx_match); % Middle of the size bin
  elseif ismember(sfield,vsize.VSDnam)
    idx_match = strcmp(vsize.VSDnam,sfield);
    info.(sfield).size_bin_st  = vsize.strnam{idx_match}; % string containing size bin limits
    info.(sfield).size_bin_mm  = vsize.sizemm(idx_match,:); % Size bin limits in mm
    info.(sfield).size_mid_mm  = vsize.sizeav_mm(idx_match); % Middle of the size bin
  end
end %% Loop through fields to allocate metadata 

%% 14 | Sort via time
[~,isort] = sort(odv.datetime);
odv = odv(isort,:);

% Also keep ctd variables if available
extras.oname = {'chloro_fluo_[mg_chl_m-3]' 'conductivity_[ms_cm-1]' 'cpar_[%]' 'depth_[m]' 'fcdom_[ppb_qse]' 'in_situ_density_anomaly_[kg_m-3]' 'nitrate_[umol_l-1]' 'oxygen_[umol_kg-1]' 'oxygen_[ml_l-1]' 'par_[umol_m-2_s-1]' 'potential_density_anomaly_[kg_m-3]' 'potential_temperature_[degc]' 'practical_salinity_[psu]' 'pressure_[db]' 'qc_flag' 'spar_[umol_m-2_s-1]' 'temperature_[degc]'};
extras.descs = {'chloro_fluo' 'conductivity' 'cpar' 'depth' 'fcdom' 'inSituDensityAnomaly' 'nitrate' 'oxygen_umol' 'oxygen_ml' 'par' 'potentialDensityAnomaly' 'potentialTemperature' 'salinity' 'pressure' 'qcFlag' 'spar' 'temperature'};
extras.units = {'mg chl/m^3'  'ms/cm'        '%'    'm'    'ppb qse' 'kg/m^3'              'umol/l'  'umol/kg'     'ml/l'     'umol/m^2/s' 'kg/m^3'            'degC'                 'psu'      'dbar'     'none'   'umol/m^2/s' 'degC'};


info =  struct();
info.size_bins = vsize;
try
  for nfield = 1:numel(odv.Properties.VariableNames)
    sfield = odv.Properties.VariableNames{nfield};
    if ismember(cols(nfield),extras.oname)
      idx_match = strcmp(extras.oname,cols(nfield));
      sfield = extras.descs{idx_match};
      odv.Properties.VariableNames{nfield} = sfield; % Replace with simple name
      info.(sfield).name      = cols{nfield};
      info.(sfield).unit      = extras.units{idx_match};
      info.(sfield).long_name = extras.descs{idx_match};
    else
      info.(sfield).name      = cols{nfield};
      info.(sfield).unit      = units{nfield};
      info.(sfield).long_name = descs{nfield};
    end
    if ismember(sfield,vsize.NSDnam)
      idx_match = strcmp(vsize.NSDnam,sfield);
      info.(sfield).size_bin_st  = vsize.strnam{idx_match};    % string containing size bin limits
      info.(sfield).size_bin_mm  = vsize.sizemm(idx_match,:);  % Size bin limits in mm
      info.(sfield).size_mid_mm  = vsize.sizeav_mm(idx_match); % Middle of the size bin
    elseif ismember(sfield,vsize.VSDnam)
      idx_match = strcmp(vsize.VSDnam,sfield);
      info.(sfield).size_bin_st  = vsize.strnam{idx_match};    % string containing size bin limits
      info.(sfield).size_bin_mm  = vsize.sizemm(idx_match,:);  % Size bin limits in mm
      info.(sfield).size_mid_mm  = vsize.sizeav_mm(idx_match); % Middle of the size bin
    end
  end %% Loop through fields to allocate metadata
catch
  fprintf('Error when trying to add details about "%s"\n',sfield)
  fprintf('Entering keyboard mode... investigate what happened\n')
  keyboard
end
end %% MAIN FUNCTION