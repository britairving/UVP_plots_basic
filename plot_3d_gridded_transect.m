function current_figure = plot_3d_gridded_transect(LON,LAT,Z,DAT,opt,bathy)
%% FUNCTION PLOT_3D_GRIDDED_TRANSECT
%  Syntax:
%    SAVE_FIGNAME = PLOT_3D_GRIDDED_TRANSECT(X,Z,DAT,opt)
%  
%  Description:
%    Creates standard contourf plot of data gridded to X and Z.
%
%  See also:
%    contourf
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Debug mode
dbstop if error

%% 1 | Configure script and contour characteristics
opt.map_proj       = 'lambert';          % map projection. Change this based on your cruise trajectory
opt.CameraPosition = [-225 -245 -22181]; % random number based on best P16N view. Change this based on your cruise trajectory
if ~isfield(opt,'data_clims')            % data color limits for colorbar and color scaling
  percentiles =  prctile(DAT(:),[0.5 99.5]);
  opt.data_clims    = percentiles;
  opt.data_clims(1) = floor(percentiles(1)/0.25)*0.25; % floor(percentfiles(1));
  opt.data_clims(2) = ceil(percentiles(2)/0.25) *0.25; % ceil(percentfiles(2));
end

%% 2 | Initialize figure
makefig; 
cmap = lansey;  % Set colormap
pause(0.2);     % pause to matlab a chance to catch up.
fig = gcf; fig.Name = opt.plot_title;
ax = gca;
hold(ax,'on'); grid(ax,'on');
title(ax, {strrep(opt.project,'_',' '); strrep(opt.plot_title,'_',' ')})
colormap(ax,cmap); 

%% 3 | Extract and format the bathymetry
% Set up projection based on available bathymetry
if exist('bathy','var')
  proj_lat = [min(min(bathy.lat)) max(max(bathy.lat))];
  proj_lon = [min(min(bathy.lon)) max(max(bathy.lon))];
  % restructure bathymetry to matrix from vector. 
  if size(bathy.lat,2) == 1 % is a vector
    fprintf('Gridding bathymetry\n')
    ylat = linspace(proj_lat(1),proj_lat(2),round(abs(diff(proj_lat))));
    xlon = linspace(proj_lon(1),proj_lon(2),round(abs(diff(proj_lon))));
    [blon,blat] = meshgrid(xlon,ylat);
    bz = griddata(bathy.lon,bathy.lat,bathy.z,blon,blat);
  else
    blat = bathy.lat;
    blon = bathy.lon;
    bz   = bathy.z;
  end
else % bathymetry was not passed through in input arguments
  llint = 1.0;
  proj_lat = [min(opt.lat)-llint max(opt.lat)+llint];
  proj_lon = [min(opt.lon)-llint max(opt.lon)+llint];
  [bz,blon,blat] = m_elev([proj_lon proj_lat]);
  bz = -bz; % + bathymetry : m_elev returns ocean depths as negative
end

%% 4 | Set map projection
m_proj(opt.map_proj,'lat',proj_lat,'lon',proj_lon,'rect','on');

%% 5 | Plot the bathymetry 
hbath = surf(ax,blon,blat,bz);
hbath.FaceColor = [0.8 0.8 0.8];
hbath.EdgeColor = 'none';
hbath.FaceAlpha = 0.8;
% lighting and shading bath...
camlight right; lighting gouraud; material dull;
hold(ax,'on');
% % highlight specific bathymetric contours
% bathy_levels   = round([0:max(opt.depth)/5:max(opt.depth)]);
% [cs,h] = contour3(ax,blon,blat,bz,bathy_levels,'color',[0.6 0.6 0.6],'LineWidth',1);
% clabel(cs,h,'fontsize',12,'color',[0.6 0.6 0.6],'LabelSpacing',1000);

%% 6 | Plot the data
try
  hscat = surf(ax,LON,LAT,Z,DAT);
  hscat.FaceLighting = 'none';
  hscat.EdgeColor = 'none';
catch
  hscat = scatter3(ax,opt.lon,opt.lat,opt.depth,50,opt.data,'filled');
end

%% 7 | Only show current section, if necessary
if isfield(opt,'section_idx')
  % grid time so can get range
  XTIME  = griddata(opt.lat,opt.depth,opt.time ,LAT,Z);
  midx   = XTIME >= opt.time(opt.section_idx(1)) & XTIME <= opt.time(opt.section_idx(end));
  % Set current visibility to transparent
  hscat.FaceAlpha = 0.2;
  % Remove data outside of section range
  LON(~midx) = nan; LAT(~midx) = nan;
  DAT(~midx) = nan; Z(~midx) = nan;
  % Replot surface plot with section data
  hscat2 = surf(ax,LON,LAT,Z,DAT);
  hscat2.FaceLighting = 'none';
  hscat2.EdgeColor = 'none';
end

%% 8 | Set up colorbar
try
  hcbar = colorbar(ax,'Location','eastoutside');
  hcbar.Label.String = strrep(opt.plot_title,'_','\_');
  % hcbar.Label.Position(1) = hcbar.Label.Position(1)*1.5; % Scootch to the right
  % hcbar.Label.Rotation = 270; % rotate text
catch
  
end
ax.CLim = opt.data_clims;      % Colorlimits
caxis(ax,opt.data_clims);      % Colorlimits

%% 9 | Update Axis characteristics
ax.ZLim = [0 ceil(max(max(bz))/100)*100]; % Surface to just below maximum bathymetry
set(ax,'ZDir','reverse','box', 'on');
% Format longitude and latitude labels
xt = xticklabels(ax);
yt = yticklabels(ax);
xt2 = format_position_labels(xt,'lon');
yt2 = format_position_labels(yt,'lat');
set(ax,'xticklabels',xt2,'yticklabels',yt2)
% set axis labels and title
ylabel(ax,'Latitude', 'fontsize', 18);
xlabel(ax,'Longitude', 'fontsize', 18);
zlabel(ax,'Depth [m]', 'fontsize', 18);
ax.CameraPosition = opt.CameraPosition;

%% 10 | Plot the coastline
% zlims = ax.ZLim;
% xlims = ax.XLim;
% ylims = ax.YLim;
% try
%   coastfile = [opt.project '_coast.mat'];
%   if ~exist(coastfile,'file')
%     m_usercoast(coastfile,'Color','k'); % to draw coast on map
%   else
%     m_coast('patch','k');
%   end
% catch
%   fprintf('could not plot coast\n')
% end
% % reset to original limits
% ax.ZLim = zlims;
% ax.XLim = xlims;
% ax.YLim = ylims;

%% 11 | Save figure
save_figname = [opt.project '_3D_transect_' opt.plot_title];
if isfield(opt,'savefig') && opt.savefig == 1
  standard_printfig( save_figname);
end

%% 12 | Return figure handle
current_figure = gcf;
try fig.Name = save_figname; end
end %% MAIN FUNCTION: 
