function plot_data_vs_depth(opt)
%% FUNCTION PLOT_DATA_VS_DEPTH
%  Syntax:
%    SAVE_FIGNAME = PLOT_DATA_VS_DEPTH(DEPTH,DATA,OPT)
%  
%  Description:
%    Waterfall plot of data that is intended to highlight vertical
%    evolution of each profile. 
%
%  See also:
%    plot_gridded_data_vs_depth_waterfall-similar plot but for gridded data
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Debug mode
% close all
dbstop if error
%% Only show current section, if necessary
if isfield(opt,'section_idx')
  depth = opt.depth(opt.section_idx);
  data  = opt.data(opt.section_idx);
else
  depth = opt.depth;
  data  = opt.data;
end
%% Defaults
cfg.split_watercolumn = 1;
cfg.maxdepth          = ceil(max(depth/50))*50; % Round maximum depth to next 50m
%% Initialize figure
makefig;
pause(0.2); % pause to matlab a chance to catch up.
fig = gcf;
fig.Name = opt.plot_title;
ax = gca;

% Find appropriate interval to separate profiles
data_percentiles = prctile(data,[1 99]);
%% Get unique profiles
profile_labels = unique(opt.profile);

%% Loop through profiles and plot data
for np = 1:numel(profile_labels) % Loop
  idx_profile = strcmp(opt.profile,profile_labels{np});
  % Plot profile data
  plot_against_depth(ax,data(idx_profile),depth(idx_profile))
  plot_against_depth(ax,smooth(data(idx_profile),3,'moving'),depth(idx_profile),ax.YLim,'g')
  ax.XLim = data_percentiles;
  title(ax,['Profile: ' strrep(profile_labels{np},'_','\_') ' | ' strrep(opt.plot_title,'_',' ')])
  if opt.savefig
    standard_printfig([profile_labels{np} '_' strrep(opt.plot_title,'_',' ')]);
  end
  pause(0.1);
  if contains(opt.plot_title,'biovolume')
    keyboard
  end
  cla(ax);
end

%% Save figure
if opt.savefig
  % Generate filename to save figure to
  save_figname = [opt.project '_profile_evolution_of_' opt.opt.plot_title];
  standard_printfig(save_figname);
end

end %% MAIN FUNCTION
