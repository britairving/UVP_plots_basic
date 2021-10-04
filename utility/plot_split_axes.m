function a_top = plot_split_axes(ax)
%% FUNCTION AX2 = PLOT_SPLIT_AXES(AX)
% Description:
%
% Syntax: AX2 = PLOT_SPLIT_AXES(AX)
%
% Input:
%
% Output:
%
% Author: Brita Irving <bkirving@alaska.edu>
%% Set split depth
split_depth = 1000; % Default to 1000m
maxdepth = ax.YLim(2);
while maxdepth < split_depth
  split_depth = ax.YLim(2)/3;
  split_depth = floor(split_depth/100)*100; % Round to 100
end
if split_depth >= maxdepth
  split_depth = floor((maxdepth/2)/100)*100;
end
%% change BOTTOM axis limits and shrink positions
% The original axis ("ax") will become the bottom axis
xlim1 = ax.XLim;
ax.YLim        = [split_depth maxdepth]; 
ax.Position(4) = ax.Position(4)/2;
ax.YAxis.Label.Units = 'normalized'; 
ax.YAxis.Label.Position(2) = 1;
ax.GridAlpha = 0.8;

%% Create an exact copy of the axis 
a_top = copyobj(ax,gcf);

%% change TOP axis limits and shrink positions
a_top.YLim = [0 split_depth]; 
if split_depth == 1000
  a_top.YTick  = [0,200,400,600,800];
end
a_top.Position(2) = a_top.Position(2) + a_top.Position(4);
a_top.YLabel      = [];
a_top.XTickLabel  = [];
a_top.XLabel      = [];
a_top.GridAlpha   = 0.8;
% make sure xlimis consistent
ax.XLim    = xlim1;
a_top.XLim = xlim1;

% resize colorbar height
if ~isempty(ax.Colorbar)
  ax.Colorbar.Position(4) = ax.Colorbar.Position(4)*2;
end

% Readjust axis so positions are correct
pause(0.1); % let plot catch up
ax.Position(3) = a_top.Position(3);
pause(0.1); % let plot catch up

end