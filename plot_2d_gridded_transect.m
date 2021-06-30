function current_figure = plot_2d_gridded_transect(X,Y,DAT,opt)
%% FUNCTION PLOT_2D_GRIDDED_TRANSECT
%  Syntax:
%    SAVE_FIGNAME = PLOT_2D_GRIDDED_TRANSECT(X,Y,DAT,opt)
%  
%  Description:
%    Creates standard contourf plot of data gridded to X and Y.
%
%  See also:
%    contourf
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Debug mode
dbstop if error


%% 1 | Configure script and contour characteristics
opt.split_watercolumn = 1;   % 1 = Highlight upper part of watercolumn, 0 = does not
opt.show_profile_loc  = 0;   % 1 = shows location (at original resolution) of each profile, 0 = does not
opt.contourlevs       = 900; % number of contour levels for plotting, higher takes more time but more detailed plot
opt.maxdepth          = ceil(max(max(Y))/50)*50; % Round maximum depth to next 50m
opt.linecolor      = [0 0 0]; % LineColor for contourf 
opt.linestyle      = 'none';  % LineStyle for contourf
if ~isfield(opt,'data_clims')
  percentiles =  prctile(DAT(:),[0.5 99.5]);
  opt.data_clims    = percentiles;
  opt.data_clims(1) = floor(percentiles(1)/0.25)*0.25; % floor(percentfiles(1));
  opt.data_clims(2) = ceil(percentiles(2)/0.25) *0.25; % ceil(percentfiles(2));
end
opt.contour_levels = [opt.data_clims(1):(diff(opt.data_clims)/opt.contourlevs):opt.data_clims(2)];

%% 2 | Initialize figure
makefig; 
cmap = lansey;  % Set colormap
pause(0.2);     % pause to matlab a chance to catch up.
fig = gcf; fig.Name = opt.plot_title;
ax = gca;
hold(ax,'on'); grid(ax,'on');
title(ax, {strrep(opt.project,'_',' '); strrep(opt.plot_title,'_',' ')})
colormap(ax,cmap); 

%% 3 | Plot the data!
[~,hc] = contourf(ax,X,Y,DAT,opt.contour_levels,...
             'LineStyle',opt.linestyle,...
             'LineColor',opt.linecolor,...
             'Fill','on');
% Adjust axis           
set(ax,'YDir','reverse','fontsize',18,'Layer','top','FontWeight','bold')
ylabel(ax,'Depth (m)')

%% 4 | Only show current section, if necessary
if isfield(opt,'section_idx')
  xmax = opt.(opt.grid_type)(opt.section_idx(end));
  xlims = [min(X(:)) xmax];
else
  xlims = [min(X(:)) max(X(:))];
end

%% 5 | Format xticklabels
ax.XLim = xlims;
xt = xticklabels(ax);
switch opt.grid_type
  case 'time'
    datetick(ax,'x','yyyy/mm/dd','KeepLimits','KeepTicks');
    xt2 = xticklabels(ax);
  otherwise 
    xt2 = format_position_labels(xt,opt.grid_type);
end
set(ax,'xticklabels',xt2,'XTickLabelRotation',45)

%% 6 | Set up colorbar
try
  hcbar = colorbar(ax);
catch
end
ax.CLim = opt.data_clims;      % Colorlimits
hcbar.Limits = opt.data_clims; % Colorlimits
caxis(ax,opt.data_clims);      % Colorlimits
hcbar.Label.String = strrep(opt.plot_title,'_',' ');

%% 7 | Mask bathymetry
try
  if min(opt.bathy) < ax.YLim(2)
    % store original depth of figure
    zbottom = ax.YLim(2);
    % Format for patch - e.g. vector into polygon
    if max(opt.bathy) > ax.YLim(2)
      bathy_basevalue = max(opt.bathy)+200;
    else
      bathy_basevalue = ax.YLim(2);
    end
    % Sort along x in case order is wrong 
    % For example, if a test profile is listed at the end, even though it's
    % time/position was at the start of the cruise)
    [sorted_x,iu] = unique(opt.(opt.grid_type));
    sorted_bathy  = opt.bathy(iu);
    
    % smooth bathymetry
    try sorted_bathy = smooth(sorted_bathy,3); end
    
    % create a closed polygon of the xaxis boundaries and bathymetry
    bathy = [bathy_basevalue; sorted_bathy; bathy_basevalue];
    xmark = [ax.XLim(1); sorted_x;  ax.XLim(end)];
    % h = plot(xmark, bathy, 'k-*','LineWidth',3); % used for testing
    fill(xmark, bathy, [0.5 0.5 0.5]);
    % reset figure maximum depth
    ax.YLim(2) = zbottom;
  end
catch
  fprintf('Could not plot bathymetry\n')
end

%% 8 | Plot locations and depth of original profiles
if opt.show_profile_loc
  dwn_sample = 1:5:numel(opt.depth);
  h = plot(ax,opt.(opt.grid_type)(dwn_sample),opt.depth((dwn_sample)),'.','Color','k','LineWidth',0.2);
end

%% 9 | highlight top of the watercolumn
if opt.split_watercolumn
  try
    % Pass output ax_top back so can change characteristics as necessary
    ax_top = plot_split_axes(ax);
    ax.Title.String = []; % Remove original title 
  end
end

%% 10 | Save figure
save_figname = [opt.project '_2D_transect_' opt.plot_title];
if isfield(opt,'savefig') && opt.savefig == 1
  standard_printfig( save_figname);
end

%% 11 | Return figure handle
current_figure = gcf;
try fig.Name = save_figname; end
end %% MAIN FUNCTION: opt

