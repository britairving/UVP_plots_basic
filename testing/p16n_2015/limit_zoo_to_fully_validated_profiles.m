function zoo = limit_zoo_to_fully_validated_profiles(validated_list,zoo)
%FUNCTION LIMIT_ZOO_TO_FULLY_VALIDATED_PROFILES
%
%  Syntax:
%    zoo = limit_zoo_to_fully_validated_profiles(validated_list,zoo)
%
%  Inputs:
%   zoo            | output from UVP_read_odv_ecotaxa_exported_zoo.m
%   validated_list | file name containing list of fully validated
%                    profiles. profile names expected to exactly match
%                    those read in from the ecotaxa file with
%                    UVP_read_odv_ecotaxa_exported_zoo.m. Simple list
%                    separated by commas. 
%                    E.g. p16n_001, p16n_006, p16n_011, p16n_012, p16n_188
%
%  Outputs:
%   zoo            | same format as input zoo, but contains only zoo data
%                    from fully validated profiles. 
% 
%  Description:
%    Reads list of fully validated profiles from a text file 
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 
fprintf('TO DO: add option to just include list of validated stations!!\n')

fprintf('------------------------\n')
fprintf('Reading fully validated profiles from: %s\n',validated_list)
try
  % Read validated profiles from file
  fileID = fopen(validated_list,'r');
  profiles = fgetl(fileID);
  profiles = strtrim(strsplit(profiles,','));
  fclose(fileID);
  % get index of matching profiles
  idx_validated = ismember(zoo.profile,profiles);
  size_profiles = size(zoo.profile,1);
  fields = fieldnames(zoo);
  for nf = 1:numel(fields)
    field = fields{nf};
    if isequal(size(zoo.(field),1),size_profiles)
      zoo.(field) = zoo.(field)(idx_validated,:);
    end
  end
catch
  fprintf('File with fully validated stations exists, but could not read profiles\n')
end