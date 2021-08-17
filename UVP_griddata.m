function [XX,YY,DAT] = UVP_griddata(x,y,d)
%FUNCTION UVP_GRIDDATA
%
%  Syntax:
%    [XX,YY,DAT] = UVP_griddata(x,y,d)
%
%  Description:
%    Grid input Nx1 data into NxM data matricies using Matlab's griddata.m
%    function. 
%
%  See also:
%    meshgrid
%    griddata
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Set up gridding characteristics
cfg.stophere = 1;         % 1 = stops at the end of the script and plots the gridded data

cfg.y_step   = 100;       % [meters] UVP data output in 5m bins from Ecotaxa
cfg.y_start  = 0;         % [meters] depth of first vertical bin
cfg.method   = 'natural'; % interpolation method 'natural','linear' or 'nearest'
cfg.maskXbyZ = 1;         % 1 = masks interpolation over vast horizontal distances, 0 = no mask applied - will interpolate over everything
cfg.scalexy  = 1;         % 1 = rescale axes so the x and y dimensions span a comparable range of values (useful for fixing horizontal streaking issues)
cfg.finite   = 1;         % 1 = ignores all NaNs - may help solve poor gridding coverage
if cfg.finite
  idx_finite = isfinite(d) & d ~= 0;
  x = x(idx_finite);
  y = y(idx_finite);
  d = d(idx_finite);
end

cfg.x_step   = abs(min(x)-max(x))/80; % degrees or days, depending on what x is
% Set the extrapolation method for scatteredInterpolant function
if cfg.maskXbyZ % set to interp method if going to mask by Z
  cfg.extrap_method = 'interp';
else% otherwise, do not extrapolate
  cfg.extrap_method = 'none';
end
%% Meshgrid X and Y
X = min(x):cfg.x_step:max(x); % x could be time, latitude, or longitude
Y = cfg.y_start:cfg.y_step:max(y);      % y is depth
[XX,YY] = meshgrid(X,Y);      % Creates NxM x and y for gridding
if cfg.scalexy 
  % Andrew McDonnells's method to solve horizontal streaking in gridding
  % Horizontal streaking arises because griddata looks for the nearest
  % neighbor for extrapolating a surface, but it does it assuming latitude
  % and depth have the same units or scaling.  As a result, it thinks a
  % point 20 degrees away at the same depth is closer than a point half a
  % degree away and only 25m difference in depth, so the gridding gets
  % stretched out in the horizontal direction.

  % Solution (see code): rescale axes so the x and y dimensions span a
  % comparable range of values, then run griddata, then convert back to
  % original units for contouring. 
  
  xscale = max(x) - min(x); % originally developed for latitude
  zscale = max(Y) - min(Y); % Depth
  xyscalefactor = xscale / zscale;
  [Xscaled,Yscaled] = meshgrid(X, Y*xyscalefactor);
end

%% Find depth boundary so griddata does not interpolate over huge areas
if cfg.maskXbyZ
  % fprintf(' Finding cruise depth boundary so can mask erroneous data - change n_interval here if masking looks strange\n')
  % within n_interval find the maximum
  n_interval = cfg.x_step; % degrees or days, depending on what x is ... may need to adjust depending on cruise length, track, etc
  data_dep   = nan(size(X));
  % Initialize
  cfg.grid_ignore = ones(size(XX));
  % Loop through each horizontal (x) location
  for jj = 1:numel(X)
    inds = find(x > (X(jj)-n_interval) & x <= (X(jj)+n_interval));
    if isempty(inds)
      continue
    end
    data_dep(jj) = max(y(inds),'omitnan');
    cfg.grid_ignore(:,jj) = YY(:,jj) > data_dep(jj);
  end
  % Convert to logical
  cfg.grid_ignore = logical(cfg.grid_ignore);
end

%% Grid data
% Try using scatteredInterpolant first
% Then use griddata if scatteredInterpolant doesn't work.
% <https://blogs.mathworks.com/loren/2009/07/22/computational-geometry-in-matlab-r2009a-part-ii/>
% <https://www.mathworks.com/matlabcentral/answers/102494-why-does-griddata-give-different-results-on-the-same-dataset>
if cfg.scalexy == 0 % Do not scale vertical distance 
  try  % scatteredInterpolant first
    ScatInterpolant = scatteredInterpolant(x,y,d,cfg.method,cfg.extrap_method);
    DAT = ScatInterpolant(XX,YY);
    cfg.gridfun = 'scatteredInterpolant';
  catch% griddata second
    [~,~,DAT] = griddata(x,y,d,XX,YY,cfg.method);
    cfg.gridfun = 'griddata';
  end
else % Scale y so values similar to x (method to reduce horizontal streaking)
  try  % scatteredInterpolant first
    ScatInterpolant = scatteredInterpolant(x,y*xyscalefactor,d,cfg.method,cfg.extrap_method);
    DAT = ScatInterpolant(Xscaled,Yscaled);
    cfg.gridfun = 'scatteredInterpolant';
  catch% griddata second
    [~,~,DAT] = griddata(x,y*xyscalefactor,d,Xscaled,Yscaled,cfg.method);
    cfg.gridfun = 'griddata';
  end
end

% Set data below depth boundary to NaNs 
if cfg.maskXbyZ
  DAT(cfg.grid_ignore) = NaN;
end

%% Pull out only finite data
% Get rid of any vertical bins containing all NaNs
% Sometimes, if there is one or two very deep profiles there will be no
% data at those depths.
idx_good_depths = find( any( isfinite(DAT), 2) );
DAT = DAT(idx_good_depths,:);
XX  = XX(idx_good_depths,:);
YY  = YY(idx_good_depths,:);

%% Print to screen the grid settings
fprintf('\n-----------------------------------------\n');
fprintf('---------UVP_griddata Settings-----------\n');
fprintf('Gridding function: %s\n',cfg.gridfun);
fprintf(' Vertical bin (m): %d\n',cfg.y_step);
fprintf(' Start depth  (m): %d\n',cfg.y_start);
fprintf(' Horizontal bin  : %f\n',cfg.x_step);
fprintf(' X/Y scaling used: %d\n',cfg.scalexy);
fprintf(' XY mask by depth: %d\n',cfg.maskXbyZ);
fprintf('-----------------------------------------\n');

%% Stop here to look at how gridding was
% Find out how good the gridding method was for this data...
percent_finite = sum(isfinite(DAT(:)))/numel(DAT)*100; % [%] of data that are not NaN
if cfg.stophere || percent_finite < 20
  makefig;ax = gca;
  pcolor(XX,YY,DAT)
  
  %contourf(XX,YY,DAT,1000,'edgecolor','none')
  shading(ax,'flat'); % Options 'interp','faceted','flat'
  set(ax,'YDir','reverse');
  hold(ax,'on');grid(ax,'on');
  axes(ax);  

  
  % show settings on plot
  text(ax,0.01,0.30,['Gridding function: ' cfg.gridfun],          'units','normalized','fontweight','bold');
  text(ax,0.01,0.25,[' Vertical bin (m): ' num2str(cfg.y_step)],  'units','normalized','fontweight','bold');
  text(ax,0.01,0.20,[' Start depth  (m): ' num2str(cfg.y_start)], 'units','normalized','fontweight','bold');
  text(ax,0.01,0.15,[' Horizontal bin  : ' num2str(cfg.x_step)],  'units','normalized','fontweight','bold');
  text(ax,0.01,0.10,[' X/Y scaling used: ' num2str(cfg.scalexy)], 'units','normalized','fontweight','bold');
  text(ax,0.01,0.05,[' XY mask by depth: ' num2str(cfg.maskXbyZ)],'units','normalized','fontweight','bold');
  
    
  title_str = {[cfg.gridfun '(' cfg.method ')  with Depth bin = ' num2str(cfg.y_step) 'm']};
  title(ax,title_str)
  ax.YLim(1) = 0;
  % Try and set colorbar
  clims = prctile(d,[0 95]);
  try
    cb = colorbar(ax); cb.Limits = clims;
    caxis(ax,cb.Limits);
  catch
    % doesn't work on my mac occasionally... it's a mystery
  end

  %% Look at the actual data with scatter
  scatter(x,y,25,d,'filled') % shows location of actual data
  ax.CLim = clims;
  % save_figname = 'p16n_griddata_example_scatteredInterpolant';
  % standard_printfig(save_figname);
end

end %% MAIN FUNCTION UVP_griddata