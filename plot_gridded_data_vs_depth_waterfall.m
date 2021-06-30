function save_figname = plot_gridded_data_vs_depth_waterfall(X,Y,DAT,opt)
%% FUNCTION PLOT_GRIDDED_DATA_VS_DEPTH_WATERFALL
%  Syntax:
%    SAVE_FIGNAME = PLOT_GRIDDED_DATA_VS_DEPTH_WATERFALL(X,Y,DAT,OPT)
%  
%  Description:
%    Waterfall plot of gridded data that is intended to highlight vertical
%    evolution of data. 
%
%  See also:
%    plot_data_vs_depth_waterfall - similar plot but for ungridded data
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Debug mode
% close all
dbstop if error

%% Defaults
cfg.split_watercolumn = 1;
cfg.maxdepth          = ceil(max(max(Y/50)))*50; % Round maximum depth to next 50m
%% Initialize figure
makefig;
pause(0.2); % pause to matlab a chance to catch up.
fig = gcf;
fig.Name = opt.plot_title;
ax = gca;
hold(ax,'on');grid(ax,'on');
title(ax,['Data evolution | ' strrep(opt.plot_title,'_',' ')])
xlabel(ax,strrep(opt.plot_title,'_',' '))
ylabel(ax,'Depth [m]')

% Find appropriate interval to separate profiles
try
  data_int = quantile(DAT(:),3); % DAT(:) because must be array, not matrix
  data_int = round(data_int(2)/2,2); % pull out middle quantile, round to 2 decimals
catch
  keyboard
end
%% Get X AXIS labels
% Y is depth
% X is latitude, longitude, or time/profiles
xunique = unique(X);
switch opt.grid_type
  case 'time'
    profile_labels = cellstr(datestr(xunique,'yyyy/mm/dd'));
    cblabel = 'Date';
  case 'lat' 
    profile_labels = strcat(cellstr(num2str(round(xunique,2))),' \circ');
    cblabel = 'Latitude';
  case 'lon'
    profile_labels = strcat(cellstr(num2str(round(xunique,2))),' \circ');
    cblabel = 'Longitude';
  otherwise
    profile_labels = cellstr(num2str(xunique,'%.2f'));
    cblabel = '';
end

% Colors of vertical profiles
clrs = lansey(numel(xunique));

%% 
depths = []; % Initialize depths array so can adjust y-axis to best visualize data
maxdat = []; % Initilaize maxdat array so can adjust x-axis to best visualize data
for np = 1:numel(xunique) % Loop
  data  = DAT(:,np);
  depth = Y(:,np); % since this is gridded, this will always be the same. 
  depths = [depths; max(depth(isfinite(data)))];
  maxdat = [maxdat; max(data+(data_int*np))];
  plot(ax,data+(data_int*np),depth,'-','LineWidth',2,'Color',clrs(np,:),'DisplayName',profile_labels{np})
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
%text(ax,0.1,0.1,[field ' offset by ' num2str(data_int) ' per profile'],'units','normalized','HorizontalAlignment','right','BackgroundColor','w')
text(ax,0.03,0.1,[strrep(opt.plot_title,'_','\_') ' offset by ' num2str(data_int) ' per profile'],'units','normalized','BackgroundColor','w')

%% Save figure
if opt.savefig
  % Generate filename to save figure to
  save_figname = [opt.project '_data_evolution_of_' opt.plot_title];
  standard_printfig(save_figname);
end
end %% MAIN FUNCTION
