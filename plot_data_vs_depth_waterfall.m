function save_figname = plot_data_vs_depth_waterfall(opt)
%% FUNCTION PLOT_DATA_VS_DEPTH_WATERFALL
%  Syntax:
%    SAVE_FIGNAME = PLOT_DATA_VS_DEPTH_WATERFALL(DEPTH,DATA,PROFILENAMES,PLOT_TITLE)
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
  profs = opt.profile(opt.section_idx);
else
  depth = opt.depth;
  data  = opt.data;
  profs = opt.profile;
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
hold(ax,'on');grid(ax,'on');
title(ax,['Profile evolution | ' strrep(opt.plot_title,'_',' ')])
xlabel(ax,strrep(opt.plot_title,'_',' '))
ylabel(ax,'Depth [m]')

% Find appropriate interval to separate profiles
data_int = quantile(data,3); % DAT(:) because must be array, not matrix
data_int = round(data_int(2)/2,4); % pull out middle quantile, round to 2 decimals
if data_int == 0
  data_int = quantile(data,3);
  data_int = data_int(2);
end
%% Get unique profiles
profile_labels = unique(profs);
% Colors of vertical profiles
clrs = jet(numel(profile_labels));

%% Initialize arrays for best visualization limits
depths = []; % Initialize depths array so can adjust y-axis to best visualize data
maxdat = []; % Initilaize maxdat array so can adjust x-axis to best visualize data
%% Loop through profiles and plot data
for np = 1:numel(profile_labels) % Loop
  idx_profile = strcmp(profs,profile_labels{np});
  % pull out max depth and maximum data for this profile
  depths = [depths; max(depth(idx_profile))];
  maxdat = [maxdat; max(data(idx_profile)+(data_int*np))];
  % Plot profile data
  plot(ax,data(idx_profile)+(data_int*np),depth(idx_profile),'-','LineWidth',2,'Color',clrs(np,:),'DisplayName',profile_labels{np})

end
ax.YLabel.String = 'Depth [m]';
ax.YDir = 'reverse';  ax.FontWeight = 'bold';

% Get an appropriate maximum depth
try
  prcntile_depths = prctile(depths,[1 90]);
  prcntile_maxdat = prctile(maxdat,[1 99]);
  ax.YLim = [0 prcntile_depths(2)];
  ax.XLim(2) = prcntile_maxdat(2);
catch
  ax.YLim = [0 cfg.maxdepth];
end

if contains(opt.plot_title,'biovolume')
  keyboard
end
% highlight top of the watercolumn
if cfg.split_watercolumn
  % Pass output ax_top back so can change characteristics as necessary
  ax_top = plot_split_axes(ax);
  % show legend
  hl = legend(ax_top,'show');
  hl.FontSize = 10;
  hl.Location = 'eastoutside';
  % Readjust axis so positions are correct
  pause(0.1); % let plot catch up
  ax.Position(3) = ax_top.Position(3);
else
  %% show legend
  hl = legend(ax,'show');
  hl.FontSize = 10;
  hl.Location = 'eastoutside';
end

% set axis labels to zero since not quantitative
ax.XTickLabel = [];
text(ax,0.03,0.1,[strrep(opt.plot_title,'_','\_') ' offset by ' num2str(data_int) ' per profile'],'units','normalized','BackgroundColor','w')


%% Save figure
if opt.savefig
  % Generate filename to save figure to
  save_figname = [opt.project '_profile_evolution_of_' opt.plot_title];
  standard_printfig(save_figname);
end
end %% MAIN FUNCTION
