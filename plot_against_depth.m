function h = plot_against_depth(ax,x,y,ylim,clr)
%% FUNCTION plot_against_depth
%
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%
if nargin < 5
  clr = 'k'; 
end
axes(ax); hold(ax,'on'); grid(ax,'on')% make sure proper axes is current
h = plot(ax,x,y,'-','color',clr,'LineWidth',2);
set(ax,'XAxisLocation','bottom','YDir','reverse');
ylabel(ax,'Depth [m]')
if nargin > 3
  set(ax,'YLim',ylim)
end

end %% FUNCTION plot_against_pressure