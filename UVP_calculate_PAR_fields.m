function [par,par_info] = UVP_calculate_PAR_fields(par,par_info)
%UVP_CALCULATE_PAR_FIELDS
%
%  Syntax:
%    UVP_calculate_PAR_fields(par,info)
%  
%  Description:
%    Calculates common variables from the standard UVP particulate output
%    from Ecotaxa. I.e. particulate abundance [#/L] and particulate
%    biovolume [ppm] or [mm^3/L].
%   
%    Calculates: 
%     datenum | Matlab datenum [Days since January 0, 0000]
%     mips    | Abundance (102-203µm)[#/L]
%     maps    | Abundance (203µm-2.58mm)[#/L]
%     tot_par_abundance  | Total Abundance (smallest-2.58mm)[#/L]
%     tot_par_abundance2 | Total Abundance (smallest->26mm)[#/L]
%     tot_par_biovolume  | Total Biovolume (smallest-2.58mm)[mm^3/L]
%     tot_par_biovolume2 | Total Biovolume (smallest->26mm)[mm^3/L]
%     NSD     | Number size distribution
%
%  Examples:
%    fedit myfun
%    fedit myfun 'input1, input2' 'output1, output2'
%
%  See also:
%    SEDIT
%    EDIT
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 1 | Create cell array with fields names (i.e. columns names) that include size limits
% Necessary because simplified fieldnames in UVP_read_[..]_par.m script
fields_wsizes = {};
fields = fieldnames(par_info);
if strcmp(fields{1},'size_bins')
  fields(1) = [];
end
for nf = 1:numel(fields)
  if strcmp(fields{nf},'size_bins')
    continue
  end
  fields_wsizes{nf} = par_info.(fields{nf}).name;
end

%% 2 | Add MATLAB's datenum
par.datenum = datenum(par.datetime);
% Add metadata
par_info.datenum.name      = 'Date';
par_info.datenum.unit      = 'Days since January 0, 0000';
par_info.datenum.long_name = 'Matlab datenum';

%% 3 | Calculate MIPS and MAPS
try
  % default mips and maps definition
  mips = [102 203];  % 'LPM_(102-128_µm)[#/L]' to 'LPM_(161-203_µm)[#/L]'
  mips_name = 'Abundance (102-203µm)[#/L]';
  %   maps = [203 26]; % 'LPM_(203-256_µm)[#/L]' to 'LPM_(>26_mm)[#/L]'
  %   maps_name = 'LPM (203µm-26mm)[#/L]';
  maps = [203 2.58]; % 'LPM_(203-256_µm)[#/L]' to 'LPM_(2.05-2.58_mm)[#/L]'
  maps_name = 'Abundance (203µm-2.58mm)[#/L]';
  idx_mips1 = find(contains(fields_wsizes,['LPM_(' num2str(mips(1)) '-']));
  idx_mips2 = find(contains(fields_wsizes,['-' num2str(mips(2)) '_µm)[#']));
  idx_maps1 = find(contains(fields_wsizes,['LPM_(' num2str(maps(1)) '-']));
  idx_maps2 = find(contains(fields_wsizes,[num2str(maps(2)) '_mm)[#']));
  if numel(idx_maps2) > 1
    idx_maps2 = idx_maps2(end);
  end
  fprintf('\n')
  fprintf(' Calculating MIPs from: %s to %s\n',fields_wsizes{idx_mips1},fields_wsizes{idx_mips2})
  fprintf(' Calculating MAPs from: %s to %s\n',fields_wsizes{idx_maps1},fields_wsizes{idx_maps2})
  % check to make sure table and column names were read correctly
  par.mips = nansum(table2array(par(:,idx_mips1:idx_mips2)),2);
  par.maps = nansum(table2array(par(:,idx_maps1:idx_maps2)),2);
  % add MIPs and MAPs to info structure 
  par_info.mips.name   = mips_name;
  par_info.maps.name   = maps_name;
  par_info.mips.unit   = '#/L';
  par_info.maps.unit   = '#/L';
  par_info.mips.long_name = ['Abundance of particles with an equivalent spherical diameter between ' num2str(mips(1)) '-' num2str(mips(2))      ' micrometers'];
  par_info.maps.long_name = ['Abundance of particles with an equivalent spherical diameter between ' num2str(maps(1)/1000) '-' num2str(maps(2)) ' millimeters'];             
catch
  fprintf(' Could NOT calculate MIPs and MAPs\n')
  keyboard
end

%% 4 | Calculate total particle abundance over the whole range of size classes <26mm
% Calculate total particle abundance from first size bin to last
% ignore fields with "=" in field (i.e. MIPS or MAPs) 
% ignore fields with ">" in field (i.e. >26_mm)
try
  for nn = 1:2
    if nn == 1 % whole range of size classes <26mm
      idx_tot_abund = (contains(fields_wsizes,'[#/L]') & ~contains(fields_wsizes,'=') & ~contains(fields_wsizes,'>'));
      name_string = 'tot_par_abundance';
    else % % whole range of size classes including those > 26mm
      idx_tot_abund = (contains(fields_wsizes,'[#/L]') & ~contains(fields_wsizes,'='));
      name_string = 'tot_par_abundance2';
    end
    
    idx_tot_1 = find(idx_tot_abund,1,'first');
    idx_tot_2 = find(idx_tot_abund,1,'last');
    % pull out string number and unit
    str_1 = erase(fields_wsizes{idx_tot_1},{'LPM_(',')[#/L]'});
    str_2 = erase(fields_wsizes{idx_tot_2},{'LPM_(',')[#/L]'});
    tmp_1 = strsplit(str_1,'-');
    num_1 = tmp_1{1};unt_1 = str_1(end-1:end);
    tmp_2 = strsplit(str_2,'-');
    if ~strcmp(str_2,'>26_mm')
      tmp_2 = strsplit(tmp_2{2},'_');
    end
    num_2 = tmp_2{1}; unt_2 = str_2(end-1:end);
    if strcmp(num_2,'>26_mm')
      num_2 = '>26';
    end
    par.(name_string) = nansum(table2array(par(:,idx_tot_abund)),2);
    tot_par_abundance_name = ['Total Particle Abundance (' num_1 unt_1 '-' num_2 unt_2 ')[#/L]'];
    fprintf(' Calculating particle abundance: %s\n',tot_par_abundance_name)
    % add total particle abundance to info structure
    par_info.(name_string).name        = tot_par_abundance_name;
    par_info.(name_string).unit        = '#/L';
    par_info.(name_string).long_name   =  ['Abundance of particles with an equivalent spherical diameter between ' num_1 unt_1 '-' num_2 unt_2];

    if strcmp(unt_1,'µm')
      num_1 = str2double(num_1)/1000;
    end
    if strcmp(num_2,'µm')
      num_2 = str2double(num_2)/1000;
    end
    par_info.(name_string).size_bin_mm = [num_1 num_2];
  end
catch
  fprintf(' Could NOT calculate total particle abundance\n')
  keyboard
end

%% 5 | Calculate total particle biovolume over the whole range of size classes <26mm
% Calculate total particle biovolume from first size bin to last
% ignore fields with "=" in field (i.e. MIPS or MAPs)
% ignore fields with ">" in field (i.e. >26_mm)
try
  for nn = 1:2
    if nn == 1 % whole range of size classes <26mm
      idx_tot_biovol = (contains(fields_wsizes,{'[ppm]' '[mm^3/L]'}) & ~contains(fields_wsizes,'=') & ~contains(fields_wsizes,'>'));
      name_string = 'tot_par_biovolume';
    else % % whole range of size classes including those > 26mm
      idx_tot_biovol = (contains(fields_wsizes,{'[ppm]' '[mm^3/L]'}) & ~contains(fields_wsizes,'='));
      name_string = 'tot_par_biovolume2';
    end
    idx_tot_1 = find(idx_tot_biovol,1,'first');
    idx_tot_2 = find(idx_tot_biovol,1,'last');
    % pull out string number and unit
    str_1 = erase(fields_wsizes{idx_tot_1},{'LPM_biovolume_(',')[ppm]'});
    str_2 = erase(fields_wsizes{idx_tot_2},{'LPM_biovolume_(',')[ppm]'});
    tmp_1 = strsplit(str_1,'-');
    num_1 = tmp_1{1};unt_1 = str_1(end-1:end);
    tmp_2 = strsplit(str_2,'-');
    if ~strcmp(str_2,'>26_mm')
      tmp_2 = strsplit(tmp_2{2},'_');
    end
    num_2 = tmp_2{1}; unt_2 = str_2(end-1:end);
    if strcmp(num_2,'>26_mm')
      num_2 = '>26';
    end
    par.(name_string) = nansum(table2array(par(:,idx_tot_biovol)),2);
    tot_par_biovolume_name = ['Total Particle Biovolume (' num_1 unt_1 '-' num_2 unt_2 ')[ppm]'];
    fprintf(' Calculating particle biovolume: %s\n',tot_par_biovolume_name)
    % add total particle biovolume to info structure
    par_info.(name_string).name        = tot_par_biovolume_name;
    par_info.(name_string).unit        = 'mm^3/L';
    par_info.(name_string).long_name   = ['Biovolume of particles with an equivalent spherical diameter between ' num_1 unt_1 '-' num_2 unt_2];
    if strcmp(unt_1,'µm')
      num_1 = str2double(num_1)/1000;
    end
    if strcmp(num_2,'µm')
      num_2 = str2double(num_2)/1000;
    end
    par_info.(name_string).size_bin_mm = [num_1 num_2];
  end 
catch
  fprintf(' Could NOT calculate total particle biovolume\n')
  keyboard
end

%% 6 | Create arrays of number size distribution & volume size distribution 
% NSD = Same as ecotaxa LPM_(size1-size2)[#/L] variables
% VSD = same as ecotaxa LPM_biovolume_(size1-size2)[ppm] variables
% ------------- Note on terminology ------------------------
% HERE - we will use VSD and NSD. 
% Different variables have been used for the same particle abundance data
% in units of number/liter, and for particle biovolume data in units of ppm
% or mm^3/liter.
% For example, NSD, CSD, PSD, or abundance  with units of[#/L]
% where NSD = number size distribution
%       CSD = concentration size distribution 
%       PSD = particle size distribution 
%       abundance = particle abundance
% For example, VSD and biovolume with units of [ppm] or [mm^3/L]
% where VSD = volume size distribution 
%       biovolume = biovolume  

nperL  = find(contains(fields_wsizes,'[#/L]')              & ~contains(fields_wsizes,'='));
biovol = find(contains(fields_wsizes,{'[ppm]' '[mm^3/L]'}) & ~contains(fields_wsizes,'='));
% Pull out NSD & VSD & convert to array
NSD = table2array(par(:,nperL));  % [#/L]
VSD = table2array(par(:,biovol)); % [mm^3/L] or [ppm]
% Biovolume fields with no value are exported from ecotaxa as 0, instead of
% blank. Therefore, replace all biovolume == 0 with NaN to avoid confusion.
VSD(VSD == 0) = NaN;

%% 7 | Remove LPM# and LPM_biovolume# variables
remove_abund_fields  = par.Properties.VariableNames(nperL);
remove_biovol_fields = par.Properties.VariableNames(biovol);
par(:,[nperL biovol])  = [];

%% 8 | Calculate DNSD (Differential number size distribution) [#/m^3/µm] & DVSD (Differential volume size distribution) [µL/m^3/µm]
% ------------- Note on terminology ------------------------
% HERE - we will use DVSD and DNSD. 
% Different variables have been used for the same data refereing to the
% noramlized particle number size distribution.
% For example, DNSD, CSDn with units of [#/L/mm] or [#/m^3/µm]
% where DNSD = differential number size distribution
%       CSDn = normalized concentration size distribution 

% Normalize size distribution to the width of the bin
fprintf(' Normalizing size distribution to the width of the bin (in milimeters) -- i.e. differential size distribution\n')
for i=1:size(NSD,2)
  wid_mm(i) = diff(par_info.size_bins.sizemm(i,:));
  DNSD(:,i) = NSD(:,i)./wid_mm(i); % [#/L/mm]                == [#/m^3/µm]
  DVSD(:,i) = VSD(:,i)./wid_mm(i); % [mm^3/L/um] or [ppm/mm] == [µL/m^3/µm]
end
% -------------------------- Notes on units --------------------------
% Convert DNSD [#/L/mm] to DNSD (Differential number size distribution) [#/m^3/µm]
% units specified by EXPORTS/SEABASS
% convert [#/L/mm] to [#/L/um] where 1mm/1000um =1
% DNSD_1 = DNSD./1000; % [#/L/um]
% convert [#/L/um] to [#/m^3/um] where 1L = 0.001m^3
% DNSD = DNSD_1./0.001; % [#/L/um]
% I.e... this cancels out because  #/L/mm is eqiuvalent to # m^-3 um^-1

% Convert VSDn [mm^3/L/mm] to DVSD (Differential volume size distribution) [µL/m^3/µm]
% 1mm^3 = 1µL
% 1L = 0.001m^3
% 1mm = 1000µm
% I.e... mm^3/L/mm is eqiuvalent to uL m^-3 um^-1
% DVSD = DVSD; % [uL m^-3 um^-1]

%% 9 | Convert to structure
par = table2struct(par,'ToScalar',true);
par.NSD = NSD;
par.VSD = VSD;
par.DNSD = DNSD;
par.DVSD = DVSD;
% Initialize tables to contain array information
par_info.NSD = table;
par_info.NSD.name = repmat({''},numel(nperL),1);
par_info.NSD.unit = repmat({''},numel(nperL),1);
par_info.NSD.long_name   = repmat({''},numel(nperL),1);
par_info.NSD.size_bin_st = repmat({''},numel(nperL),1);
par_info.NSD.size_bin_mm = nan(numel(nperL),2);
par_info.NSD.size_mid_mm = nan(numel(nperL),1);
% Initialize tables to contains the same fields, info will be overwritten
par_info.VSD = par_info.NSD;
% Loop through fields (could also use number of size bins)
for nn = 1:numel(remove_abund_fields)
  % Abundance fields
  field = remove_abund_fields{nn};
  par_info.NSD.name{nn} = strrep(par_info.(field).name,'LPM_','Abundance ');
  par_info.NSD.unit{nn} = par_info.(field).unit;
  par_info.NSD.long_name{nn}   = par_info.(field).long_name;
  par_info.NSD.size_bin_st{nn} = par_info.(field).size_bin_st;
  par_info.NSD.size_bin_mm(nn,:) = par_info.(field).size_bin_mm;
  par_info.NSD.size_mid_mm(nn)   = par_info.(field).size_mid_mm;
  % Biovolume fields
  field = remove_biovol_fields{nn};
  par_info.VSD.name{nn} = strrep(par_info.(field).name,'LPM_biovolume_','Biovolume ');
  par_info.VSD.unit{nn} = par_info.(field).unit;
  par_info.VSD.long_name{nn}     = par_info.(field).long_name;
  par_info.VSD.size_bin_st{nn}   = par_info.(field).size_bin_st;
  par_info.VSD.size_bin_mm(nn,:) = par_info.(field).size_bin_mm;
  par_info.VSD.size_mid_mm(nn)   = par_info.(field).size_mid_mm;
end
% Replicate to differential and change units, names as necessary
par_info.DNSD = par_info.NSD;
par_info.DNSD.name = strrep(par_info.DNSD.name,'[#/L]','[#/L/mm]');
par_info.DNSD.name = strcat('Normalized',{' '}, par_info.DNSD.name);
par_info.DNSD.unit = strrep(par_info.DNSD.unit,'#/L','#/L/mm');
par_info.DNSD.long_name = strcat(par_info.DNSD.long_name,' normalized to the width of the bin');
% Replicate to differential and change units, names as necessary
par_info.DVSD = par_info.VSD;
par_info.DVSD.name = strrep(par_info.DVSD.name,{'[ppm]'},'[mm^3/L/mm]');
par_info.DVSD.name = strcat('Normalized',{' '}, par_info.DVSD.name);%strrep(par_info.DVSD.name,'LPM_','nLPM_');
par_info.DVSD.unit = strrep(par_info.DVSD.unit,'#/L','#/L/mm');
par_info.DVSD.long_name = strcat(par_info.DVSD.long_name,' normalized to the width of the bin');

% Remove fields from par_info structure
par_info =  rmfield(par_info,remove_abund_fields);
par_info =  rmfield(par_info,remove_biovol_fields);

%% 10 | Calculate mean size
% Want the mean size of particles in each image. 
fprintf(' Calculating mean size\n')
% Initialize meansize
par.meansize = nan(size(par.profile));
% Loop through each row and calculate the mean size of particles in that
% row (i.e. in each profile, at each depth)
for irow = 1:numel(par.profile)
  NSD_row = par.NSD(irow,:);                         % [#/L]
  N_row   = par.NSD(irow,:)*par.SampledVolume(irow); % [#]
  % Now we have the number of particles in each size bin at this specific
  % profile depth
  % Next, create array containing an individual entry of the particle size
  % for as many particles there are in that size bin
  allsz = [];
  for nsize = 1:numel(par_info.size_bins.sizeav_mm)
    if isnan(N_row(nsize))
      continue
    end
    try
    allsz = [allsz; ones(round(N_row(nsize)),1)*par_info.size_bins.sizeav_mm(nsize)];
    catch
      keyboard
    end
  end
  par.meansize(irow) = nanmean(allsz);
end
% Add metadata
par_info.meansize.name      = 'Particle mean size [mm]';
par_info.meansize.unit      = 'mm';
par_info.meansize.long_name = 'Mean size of particles in each image';

%% 11 | Calculate slope of size distribution with a linear fit
fprintf(' Calculating slope of size distribution with a linear fit\n')
% Do not include the last size class (>26mm) because the bin center is NaN
idx_use = isfinite(par_info.size_bins.sizeav_mm);
% Initialize slobe of the size distribution
par.slope_b = nan(size(par.profile));
for irow = 1:numel(par.profile)
  x = log10(par_info.size_bins.sizeav_mm(idx_use))';
  y = log10(par.DNSD(irow,idx_use));
  % get rid of infinite and nan values, rewrote 12/22/2020
  remove_inf_nan = isinf(y) | isnan(x) | isnan(y);
  x(remove_inf_nan) = [];
  y(remove_inf_nan) = [];
  % linear regression 
  try
    pf = polyfit(x,y,1);
    if isnan(pf(1))
      fprintf(' slope of the size distribution linear fit did not work..\n')
      keyboard
    end
    par.slope_b(irow) = pf(1);
  catch
    fprintf(' slope of the size distribution linear fit did not work..\n')
    keyboard
  end
end
% Add metadata
par_info.slope_b.name      = 'Slope of Particle size distribution';
par_info.slope_b.unit      = '';
par_info.slope_b.long_name = 'Slope of size distribution calculated with a linear fit';

%% 12 | Calculate relative uncertainty due to counting statistics
fprintf(' Calculating relative uncertainty based on counting statistics\n')
% relative uncertainty due to counting statistics
% --- Discussion around uncertainty during SeaBASS submission Feb 2021 ---
% See https://www.oceanopticsbook.info/view/optical-constituents-of-the-ocean/level-3/creating-particle-size-distributions-data
% For example based, for the large bins, on counting statistics (e.g.
% \sqrt(N)) and at the small end, based also on uncertainties in size. A
% single pixel can include particles smaller than a single pixel but that
% scatter enough to register. Thus, a small particle sitting in between
% pixels can be significantly overestimated in size. This over-estimate
% (+1pixel in every direction) gets relatively smaller the larger the
% particle. Thresholding also induces uncertainties in size, as it defines
% when the signal is large enough to register.
% Since the smallest size bin we are considering is 50.8-64µm and a pixel
% is ~0.01µm, ignore the contribution to uncertainty based on size.

% Due to no published uncertainties for UVP particle abundance or biovolume
% calculations, we decided to use relative uncertainty for the time being.
par.rel_unc = sqrt( par.NSD .* par.SampledVolume) ./ (par.NSD .* par.SampledVolume); % units of percent [%]
% Add metadata
par_info.rel_unc.name      = 'Relative uncertainty';
par_info.rel_unc.unit      = '%';
par_info.rel_unc.long_name = 'relative_unc = sqrt(NSD*SampledVolume)/(NSD*SampledVolume)';

% % Calculate uncertainty in NSD
% par.u_NSD  = par.rel_unc .* par.NSD;  % [%]*[#/L]   = [#/L]       can convert to [#/m^3] later
% par_info.u_NSD.name      = 'NSD uncertainty [#/L]'; 
% par_info.u_NSD.unit      = '#/L';
% par_info.u_NSD.long_name = 'Uncertainty in the number size distribution [#/L] where NSD_unc = relative_unc*NSD';
% % Calculate uncertainty in VSD
% par.u_VSD  = par.rel_unc .* par.VSD;  % [%]*[mm^3/L] = [mm^3/L]    canconvert to [µL/m^3] later
% par_info.u_VSD.name      = 'VSD uncertainty [mm^3/L]';
% par_info.u_VSD.unit      = 'mm^3/L';
% par_info.u_VSD.long_name = 'Uncertainty in the volume size distribution [mm^3/L] or [ppm] where VSD_unc = relative_unc*VSD';
% % Calculate uncertainty in DNSD
% par.u_DNSD = par.rel_unc .* par.DNSD; % [%]*[#/L/m] = [#/L/m] = [#/m^3/µm]
% par_info.u_DNSD.name      = 'DNSD uncertainty [#/L/mm]';
% par_info.u_DNSD.unit      = '#/L/mm';
% par_info.u_DNSD.long_name = 'Uncertainty in the differential number size distribution [#/L/mm] or [#/m^3/µm] where DNSD_unc = relative_unc*DNSD';
% % Calculate uncertainty in DVSD
% par.u_DVSD = par.rel_unc .* par.DVSD; % [%]*[mm^3/L/mm] = [mm^3/L/mm] = [µL/m^3/µm]
% par_info.u_DVSD.name      = 'DVSD uncertainty [mm^3/L/mm]';
% par_info.u_DVSD.unit      = 'mm^3/L/mm';
% par_info.u_DVSD.long_name = 'Uncertainty in the differential volume size distribution [mm^3/L/mm] or [µL/m^3/µm] where DVSD_unc = relative_unc*DVSD';

end %% MAIN FUNCTION PAR = UVP_CALCULATE_PAR_FIELDS(PAR)