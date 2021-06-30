function data = uvp_remove_last_two_depth_bins(data)
%%  UVP_REMOVE_LAST_TWO_DEPTH_BINS
%
%  Syntax:
%    data = uvp_remove_last_two_depth_bins(data)
%
%  Description:
%    Removes the two deepest depth bins from each profile. 
%    Since the UVP is installed on a CTD rossette there is more time spent
%    at the deepest point in the profile, as the winch is stopped and
%    reverse. Because of this, the counts in the lowest depth bins can be
%    biased to higher counts. 
%
%  Inputs:
%    data = structure containing UVP (zoo or par) or CTD data
%
%  Outputs:
%    data = structure containing UVP (zoo or par) or CTD data with lowest
%    two depth bins removed. 
%
%  Author: Brita Irving <bkirving@alaska.edu>
%% Pull out unique profiles
profiles = unique(data.profile,'stable');
%% Loop through each profile and identify the two deepest depth bins
remove_these_depths = [];
for np = 1:numel(profiles)
  idx_p = find(strcmp(data.profile,profiles{np}));
  [~,deepest_zbin] = max(data.Depth(idx_p));
  if deepest_zbin == 1
    idx_2deepest = idx_p(deepest_zbin:deepest_zbin+1);
  else
    idx_2deepest = idx_p(deepest_zbin-1:deepest_zbin);
  end
  % remove the deepest, and next deepest zbin
  remove_these_depths = [remove_these_depths; idx_2deepest];
end
% Visually confirm this method is accurate
% figure,plot(data.latitude,data.Depth,'ko')
% hold on
% plot(data.latitude(remove_these_depths),data.Depth(remove_these_depths),'g*')

%% Loop through data fields and remove the two deepest depth bins
fields = fieldnames(data);
depth_size = size(data.Depth,1);
for nf = 1:numel(fields)
  field = fields{nf};
  if isequal(size(data.(field),1),depth_size)
    data.(field)(remove_these_depths,:) = [];
  end
end

end %% MAIN FUNCTION