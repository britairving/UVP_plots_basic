function [zoo,info] = UVP_read_odv_ecotaxa_exported_zoo(zoo_file)
%% function UVP_read_odv_ecotaxa_exported_zoo
%
%  Syntax:
%    [zoo,zoo_info] = UVP_read_odv_ecotaxa_exported_zoo(zoo_file)
%
%  Description:
%    Reads detailed ZOO ODV file exported from Ecotaxa
%
%  Inputs:
%    zoo_file = Path to detailed ODV file exported from Ecotaxa
%
%  Outputs:
%    zoo = table containing data exported from Ecotaxa
%
%  Dependencies:
%   Matlab R2016b or later (uses detectImportOptions.m and contains.m)
%
%  Notes:
%     This script uses subfunctions for organization.
%     global variables (blue) are accessible across subfunctions.
%
%  Author: Brita Irving <bkirving@alaska.edu>
%% 0 | See if par_file exists
% File must have been exported from Ecotaxa Particle Module in Detailed ODV
% format.  <https://ecotaxa.obs-vlfr.fr/part/>
% IF file does not exist, or was not passed through as an input argument
if nargin < 1
  error('Must pass through full file path to detailed ODV file exported from Ecotaxa')
elseif ~exist(zoo_file,'file')
  fprintf(' ** File not found: %s **\n',zoo_file)
  fprintf(' File must have been exported from Ecotaxa Particle Module in Detailed ODV format\n')
  error('  Could not find par file. Expected full file path to detailed ODV file exported from Ecotaxa')
end

%% 1 | Read column names and header information from ODV file
% delimiter can be ';' or '\t' (tab)
table_options = detectImportOptions(zoo_file);
table_options.CommentStyle = '//';
fprintf('\n  Reading ODV ZOO data from... %s\n',zoo_file)
% Read header
% column names usually at line 7, but go to line 20 in case extra
% header information is in file
hdr = cell(20,1);
fileID = fopen(zoo_file);
for nline = 1:20
  hdr{nline} = fgetl(fileID);
end
fclose(fileID);
column_hdrline = find(contains(hdr,'Cruise'));
% resize header
hdr = hdr(1:column_hdrline);
cols = strsplit(hdr{column_hdrline},table_options.Delimiter); % split hdr by odv delimiter
% store original column names
cols_orig = cols;
for ic = 1:numel(cols)
  cols{ic} = strrep(cols{ic},' ','_');
  cols{ic} = strrep(cols{ic},'__','_'); % clean up unnecessary underscores
end

% ecotaxa changed the format of the unit part of column names
% this happend sometime in 2019 (I think)
% old: [#/m3]     [ppm]
% new: [# m-3]    [mm3 l-1]

% For simplicity, change column names back to original format
cols = strrep(cols,'[#_m-3]','[#/m^3]');   % abundance
cols = strrep(cols,'[mm3_l-1]','[mm^3/L]');% biovolume
cols = strrep(cols,'[ppm]'    ,'[mm^3/L]');% biovolume

%% 2 | Reformat variable names to be more matlab friendly
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

%% 3 | Read ODV ZOO data
odv = readtable(zoo_file,table_options);

%% 4 | Fill meta variables
% ODV format has a single row of metadata per station. 
fnames = fieldnames(odv);
for iv = 1:numel(meta_idx)
  vname = fnames{meta_idx(iv)};
  odv.(vname) = fillmissing(odv.(vname),'previous');
end % loop through meta variables

%% 5 | Pull out indices of abundance, biovolume, and average ESD data
idx_abund = find(contains(cols,'[#/m^3]')); % abundance
idx_biovl = find(contains(cols,'[mm^3/L]'));% biovolume
idx_avesd = find(contains(cols,'[mm]'));    % avg esd
if ~isequal(numel(idx_abund),numel(idx_biovl),numel(idx_avesd))
  error('could not pull out fields correctly using the expected units')
end

%% 6 | Add metadata for each variable in a structure
descs  = {'Cruise' 'Station' 'Profile' 'Dataowner' 'UVP_Rawfilename' 'Instrument' 'SerialNumber' 'CTDrosettefilename' 'Datetime_UTC'        'Latitude'      'Longitude'   'Depth' 'Sampled_Volume'}; 
units  = {'none'   'none'    'none'    'none'      'none'            'none'       'none'         'none'               'dd-mm-yyyy HH:MM:SS' 'degrees_north' 'degrees_east' 'm'    'L'};
info =  struct();
for nfield = 1:numel(descs)
  sfield = odv.Properties.VariableNames{nfield};
  info.(sfield).name         = cols{nfield};
  info.(sfield).unit         = units{nfield};
  info.(sfield).long_name    = descs{nfield};
  info.(sfield).ecotaxa_name = cols_orig{nfield};
end %% Loop through fields to allocate metadata 

%% 7 | Parse fieldnames for taxonomic names
names.original = erase(cols(idx_abund),'_[#/m^3]');
% initialize parent and child taxonomic names
names.child  = repmat({''},size(names.original));
names.parent = repmat({''},size(names.original));
for ntaxa = 1:numel(names.original)
  str = strsplit(names.original{ntaxa},'(');
  names.child(ntaxa)  = str(1);
  % If there is a parent name 
  if numel(str) == 2
    names.parent(ntaxa) = erase(str(2),')');
  end
end

%% 8 | Simplify fieldnames
% This is not 100% necessary because these fields are moved into an array
% so fieldnames are not preserved, but still make it easier for viewing at
% this stage.
% prefix fieldname with type
odv.Properties.VariableNames(idx_abund) = strcat('abund_',odv.Properties.VariableNames(idx_abund));
odv.Properties.VariableNames(idx_biovl) = strcat('biovol_',odv.Properties.VariableNames(idx_biovl));
odv.Properties.VariableNames(idx_avesd) = strcat('avgesd_',odv.Properties.VariableNames(idx_avesd));

% remove type & unit appended to the fieldname
odv.Properties.VariableNames(idx_abund) = erase(odv.Properties.VariableNames(idx_abund),{'___M_3_' '__M_3_'});
odv.Properties.VariableNames(idx_biovl) = erase(odv.Properties.VariableNames(idx_biovl),{'_Biovolume_mm3L_1_' 'Biovolume_mm3L_1_'});
odv.Properties.VariableNames(idx_avesd) = erase(odv.Properties.VariableNames(idx_avesd),{'_Avgesd_mm_' 'Avgesd_mm_'});

%% 9 | Create tables to contain metadata about abundance, biovolume, and average ESD data
% Initialize tables to contain array information
info.abund = table;
info.abund.name         = repmat({''},numel(names.original),1);
info.abund.unit         = info.abund.name;
info.abund.long_name    = info.abund.name;
info.abund.ecotaxa_name = info.abund.name;
info.abund.child        = info.abund.name;
info.abund.parent       = info.abund.name;
% Initialize biovolume and avgesd tables too
info.biovol = info.abund;
info.avgesd = info.abund;

% Loop through taxonomic names and store associated information
for ntaxa = 1:numel(names.original)
  abund_variable =  erase(odv.Properties.VariableNames{idx_abund(ntaxa)},'abund_');% remove the abund_
  child_parent   = strsplit(abund_variable,'_');
  child          = child_parent{1};
  % abundance
  info.abund.name{ntaxa}         = [child '[#/m^3]'];
  info.abund.unit{ntaxa}         = '#/m^3';
  info.abund.long_name{ntaxa}    = char(strcat('Abundance of',{' '},names.original{ntaxa}));
  info.abund.ecotaxa_name{ntaxa} = cols_orig{idx_abund(ntaxa)};
  info.abund.child{ntaxa}        = names.child{ntaxa};
  info.abund.parent{ntaxa}       = names.parent{ntaxa};
  % biovolume
  info.biovol.name{ntaxa}         = [child '[mm^3/L]'];
  info.biovol.unit{ntaxa}         = 'mm^3/L'; % [mm^3/L] or [ppm]
  info.biovol.long_name{ntaxa}    = char(strcat('Biovolume of',{' '},names.original{ntaxa}));
  info.biovol.ecotaxa_name{ntaxa} = cols_orig{idx_biovl(ntaxa)};
  info.biovol.child{ntaxa}        = names.child{ntaxa};
  info.biovol.parent{ntaxa}       = names.parent{ntaxa};
  % average equivalent spherical diameter 
  info.avgesd.name{ntaxa}         = [child '[mm]'];
  info.avgesd.unit{ntaxa}         = 'mm'; 
  info.avgesd.long_name{ntaxa}    = char(strcat('Average equivalent spherical diameter of',{' '},names.original{ntaxa}));
  info.avgesd.ecotaxa_name{ntaxa} = cols_orig{idx_avesd(ntaxa)};
  info.avgesd.child{ntaxa}        = names.child{ntaxa};
  info.avgesd.parent{ntaxa}       = names.parent{ntaxa};
end

%% 10 | Query WoRMS to pull out AphiaID match for each taxonomic name
% World Register of Marine Species (WoRMS) is the official taxonomic
% reference list for the Ocean Biodiversity Information System (OBIS) and
% can be used to provide machine-readable taxonomic identifiers for living
% organisms.  Taxonomic names are matched with a WoRMS AphiaID to the
% lowest taxonomic rank that can be identified.
try 
  save_taxa_filename = strrep(zoo_file,'.txt','_taxa.mat');
  if ~exist(save_taxa_filename,'file')
    % Initialize table used as input to WoRMS_AphiaID_taxa_match.m script.
    % Inputs: taxa | table with fields Name & Name_parent
    taxa = table();
    taxa.Name          = names.child';  % taxa name, lowest level.        Required field in WoRMS_AphiaID_taxa_match.m.
    taxa.Name_parent   = names.parent'; % taxa parent name, if available. Required field in WoRMS_AphiaID_taxa_match.m.
    
    check_children_manual = false; % true = checks full classifcation of parent if child doesn't have a match and asks user. false = skips this.
    fprintf('  Querying WoRMS database for AphiaID\n')
    taxa = WoRMS_AphiaID_taxa_match(taxa,check_children_manual,'ecotaxa');
    fprintf('  Saving taxa matches to %s\n',[pwd filesep save_taxa_filename])
    save(save_taxa_filename,'taxa');
  else % Load data instead of rerunning everything
    fprintf('  Loading..%s\n',save_taxa_filename)
    load(save_taxa_filename);
  end
  % Store AphiaIDs to zoo information 
  info.abund.AphiaID  = repmat({''},numel(names.original),1);
  info.biovol.AphiaID = repmat({''},numel(names.original),1);
  info.avgesd.AphiaID = repmat({''},numel(names.original),1);
  info.abund.AphiaID_parent  = repmat({''},numel(names.original),1);
  info.biovol.AphiaID_parent = repmat({''},numel(names.original),1);
  info.avgesd.AphiaID_parent = repmat({''},numel(names.original),1);
  for ntaxa = 1:numel(names.original)
    info.abund.AphiaID(ntaxa)  = taxa.AphiaID(ntaxa);
    info.biovol.AphiaID(ntaxa) = taxa.AphiaID(ntaxa);
    info.avgesd.AphiaID(ntaxa) = taxa.AphiaID(ntaxa);
    info.abund.AphiaID_parent(ntaxa)  = taxa.AphiaID_parent(ntaxa);
    info.biovol.AphiaID_parent(ntaxa) = taxa.AphiaID_parent(ntaxa);
    info.avgesd.AphiaID_parent(ntaxa) = taxa.AphiaID_parent(ntaxa);
  end
catch
  fprintf(' WoRMS AphiaID query did not work\n')
end

%% 11 | Pull out arrays of data types & Remove individual abundance, biovolume, and avgesd variables
abund  = table2array(odv(:,idx_abund));
biovol = table2array(odv(:,idx_biovl));
avgesd = table2array(odv(:,idx_avesd));
% Keep all taxa varaibles in arrays for ease of use
odv(:,[idx_abund idx_biovl idx_avesd])  = [];

%% 12 | Sort via time
[~,isort] = sort(odv.datetime);
odv = odv(isort,:);

%% 13 | Convert to structure 
zoo = table2struct(odv,'ToScalar',true);

%% 14 | Add MATLAB's datenum
zoo.datenum = datenum(zoo.datetime);
% Add metadata
info.datenum.name      = 'Date';
info.datenum.unit      = 'Days since January 0, 0000';
info.datenum.long_name = 'Matlab datenum';

%% 15 Add NxM arrays of abundance, biovolume, and avgesd
zoo.abund = abund;
zoo.biovol = biovol;
zoo.avgesd = avgesd;
end %% MAIN FUNCTION