function plot_general_map(ax,lon,lat,region)
%% function plot_general_map
%
%  Syntax:
%    map_fig = plot_general_map(ax,lon,lat)
%
%  Description:
%    Generate basic map with bathymetry
%
%  Inputs:
%
%  Notes:
%     Coastline:  If coastfile does not exist, writes coastfile to pwd
%     Bathymetry: usings m_gshhs_f to generate bathymetry - edit with detailed
%                 bathymetry if desired
%
%  Author: Brita Irving <bkirving@alaska.edu>
%
%%
axes(ax);

colormap gray
dbstop if error 
% maxy = ceil(max(meta.bottomdepth)/50)*50;
% cblimits = [-maxy 0]; % bathymetry
% catch case where longitude and latitude are not capitalized

% catch if straddles 180 degrees E
if ~all(lon < 0) || ~all(lon > 0)
  lon = wrapTo360(lon);
end
la = 0.5;
lo = 0.5;
proj_lat = [min(lat)-la max(lat)+la];
proj_lon = [min(lon)-lo max(lon)+lo];

m_proj('lambert','lat',proj_lat,'lon',proj_lon,'rect','on');
% generate generic elevation and coastline with m_map
coastfile = [region '_coast.mat'];
if exist(coastfile,'file')
  m_usercoast(coastfile,'Color','k');
else
  m_gshhs_f('save',[region '_coast']);
end

hold on
m_grid('xtick', 9,'ytick', 8,'fontsize', 14);

m_plot(lon,lat,'ok','MarkerFaceColor','k','MarkerSize',6,'clipping','on');

ylabel('Latitide')
xlabel('Longitude')
title(strrep(region,'_','\_'))
% set(gca,'XDir','reverse')
set(gca,'ydir','normal')
% set colorbar for bathymetry
% cb = colorbar;
% cb.Limits = cblimits;
% cb.Label.String = 'Bathymetry [m]';
% cb.FontAngle = 'italic';
% caxis(cblimits);
% labels = -str2num(char(cb.TickLabels));
% cb.TickLabels = num2str(labels);

end %% FUNCTION plot_station_map