function plot_map_inset(fig,opt)
%% function plot_map_inset
%
%  Syntax:
%    plot_map_inset(fig,opt)
%
%  Description:
%    Plot a small map in the upper left corner of the current figure
%
%  Inputs:
%    fig | current graphical figure (i.e. gcf)
%    opt | structure containing 'lat' and 'lon'
%
%  Author: Brita Irving <bkirving@alaska.edu>
%% 0 | Check that all the information is available
if ~isfield(opt,'lat') || ~isfield(opt,'lon')
  error('Need to pass "lon" and "lat" as fields to the options structure');
end

%% 1 | Create another axes
fig.Units = 'normalized';

% Position the axis in the top left corner
amap = axes('Position',[0.01 0.75 0.1 0.22],'FontSize',10);

%% Set map projection
% Get map projection of current figure, if possible
axes(amap); 
ll_int = 1;
m_proj('lambert','lat',[min(opt.lat)-ll_int max(opt.lat)+ll_int],'lon',[min(opt.lon)-ll_int max(opt.lon)+ll_int],'rect','on');


%% Plot coast
try
  coastfile = [opt.project '_coast.mat'];
  if ~exist(coastfile,'file')
    m_usercoast(coastfile,'Color','k'); % to draw coast on map
  else
    m_gshhs_c('patch',script_opt.land_color); % Course resolution
  end
catch
  %fprintf('could not plot coast\n')
end

%% Plot cruise trajectory
axes(amap);% make sure the small axis is current
hold(amap,'on')
m_grid('box','fancy','tickdir','in','XaxisLocation','top');
m_plot(opt.lon,opt.lat,'ko','MarkerFaceColor','k','MarkerSize',3);

%% Highlight current range (i.e. if looping through sections)
if isfield(opt,'section_idx')
  m_plot(opt.lon(opt.section_idx),opt.lat(opt.section_idx),'go','MarkerFaceColor','g','MarkerSize',3);
end
%% Save figure
if isfield(opt,'savefig') && opt.savefig == 1
  save_figname = [fig.Name '_map'];
  standard_printfig( save_figname);
end

end %% FUNCTION plot_station_map