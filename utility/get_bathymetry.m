function bathy = get_bathymetry(lon,lat,project_name)
%% function get_bathymetry
%
%  Syntax:
%    bathy = get_bathymetry(lon,lat,project_name)
%
%  Description:
%    Save coastfile and loads & saves bathymetry using input longitude,
%    latitude. Generates latitude and longitude bounding boxes with lat_int
%    and lon_int defined below. 
%    Bathymetry is read in from file [project_name '_bathymetry.txt'] if
%    possible, if not, it is loaded using m_elev. 
%    See below for instructions on a simple way to download high resolution
%    bathymetry data from the web. 
% 
%  Inputs:
%    lon          | vector data with longitudes in decimal degrees. 
%    lat          | vector data with latitudes in decimal degrees. 
%    project_name | char or string with shorthand for project/dataset
%
%  Outputs:
%    bathy | structure with fields lon, lat, z, and bathy_transect
%            bathy.lon is the same format as input [0 360] or [-180 180]
%            bathy.lat is the latitude
%            bathy.z   is the bathymetry depth (+ in the ocean)
%            bathy.bathy_transect is the "z" at each input lat/lon
%
%  Notes:
%     Coastline:  If coastfile does not exist, writes coastfile to pwd
%     Bathymetry: usings m_elev to generate bathymetry 
%    
%     ***** YOU MUST HAVE M_MAP TOOLBOX INSTALLED!*****
%         <https://www.eoas.ubc.ca/~rich/map.html>
%
%  Author: Brita Irving <bkirving@alaska.edu>
%% Check that all arguments are available
if nargin < 3
  error('Must pass in longitudes, latitudes, and project_name name/shorthand')
end
%% Longitude to [0 360] format
lon_orig = lon;
lon = wrapTo360(lon);

%% Define latitude and longitude bounds
lat_int = 1.5;
lon_int = 1.5;
proj_lat = [min(lat)-lat_int max(lat)+lat_int];
proj_lon = [min(lon)-lon_int max(lon)+lon_int];
if ~isequal(lon,lon_orig)
  proj_lon = wrapTo180(proj_lon);
end
%% Generate full resoltuion coastline data with m_map
try
  coastfile = fullfile(pwd,[project_name '_coast.mat']);
  if ~exist(coastfile,'file')
    % Create projection so can generate coastline
    fprintf('Saving coastline data to file: %s\n', coastfile)
    m_proj('lambert','lat',proj_lat,'lon',proj_lon,'rect','on');
    m_gshhs_f('save',coastfile);
    %  m_usercoast(coastfile,'Color','k'); % to draw coast on map
  end
end

%% Read full resoltuion elevation with m_map
% How to download your own high resolution bathymetric data
% 1. EXTRACT XYZ GRID - TOPOGRAPHY OR GRAVITY 
%    <https://topex.ucsd.edu/cgi-bin/get_data.cgi>
% 2. COPY AND SAVE WEB QUERY RESULTS TO:  fullfile(project_folder, 'bathymetry.txt');
%    Longitude Latitude Elevation
try
  bathy_file =  fullfile(pwd, [project_name '_bathymetry.txt']);
  bathymetry = readtable(bathy_file);
  % Expect longitude, latitude, elevation (negative)
  if size(bathymetry,2) == 3
    bathy.lon = table2array(bathymetry(:,1)); % LONGITUDE
    bathy.lat = table2array(bathymetry(:,2)); % LATITUDE
    bathy.z   = table2array(bathymetry(:,3)); % BATHYMETRY 
    % change elevation (- in ocean) to bathymetry (+ in ocean)
    bathy.z  = -1.0*bathy.z ;
  else
    fprintf('Unknown format in bathymetric file, add specific case here\n')
    keyboard
  end
catch
  try
    fprintf('No high resolution bathymetry file exists, loading bathymetry with m_elev\n')
    [ELEV,LONG,LAT] = m_elev([proj_lon proj_lat]);
    % Rearrange to 
    bathy.lon = LONG(:); % Longitude in decimal degrees
    bathy.lat = LAT(:);  % Latitude in decimal degrees
    bathy.z   = ELEV(:); % Returns negative elevations in the ocean
    % change elevation (- in ocean) to bathymetry (+ in ocean)
    bathy.z  = -1.0*bathy.z ;
  catch % this will fail if user does not have m_map installed
    error('Could not load or generate bathymetry...')
    % The best way to do this is to read in high resoltuion bathymetry
    % How to download your own high resolution bathymetric data
    % 1. EXTRACT XYZ GRID - TOPOGRAPHY OR GRAVITY
    %    <https://topex.ucsd.edu/marine_topo/>
    % 2. ENTER YOUR LATITUDE AND LONGITUDE BOUNDS, AND SELECT TOPOGRAPHY
    % 3. COPY AND SAVE WEB QUERY RESULTS TO: [project_name '_bathymetry.txt']
    %    Longitude Latitude Elevation
  end
end

%% Convert longitude to [-180 180] if that was the original format
if ~isequal(lon,lon_orig)
  bathy.lon = wrapTo180(bathy.lon);
end

%% Pull out bathymetry along cruise track
try
  % Loop through each latitude and longitude point and pull out the nearest
  % bathymetry
  fprintf('Pulling out nearest bathymetry along cruise track\n')
  
  % Latitude/longitude at each profile are the same, so only look at the
  % unique entries
  [ulats,iu] = unique(lat);
  ulons = lon_orig(iu);
  bathy.bathy_transect = nan(size(ulats));
  bathy.lat_transect   = nan(size(ulats));
  bathy.lon_transect   = nan(size(ulats));

  try progressbar; end% Init single bar
  for nn = 1:numel(ulats)
    % get the closest distance
    lons = repmat(ulons(nn),size(bathy.lon));
    lats = repmat(ulats(nn),size(bathy.lon));
    % Distance between points on sphere or ellipsoid
    % Returns arclen in degrees
    [arclen,~] = distance(lats, lons, bathy.lat, bathy.lon);
    [~,idxmin] = min(arclen);
    % Pull out the closest bathymetry
    bathy.bathy_transect(nn) = bathy.z(idxmin);
    bathy.lat_transect(nn)   = ulats(nn);
    bathy.lon_transect(nn)   = ulons(nn);
    clearvars arclen idxmin lats lons
    % Update progressbar
    try progressbar(nn/numel(ulats)); end
  end
catch
  fprintf('could not pull out bathymetry along cruise track\n')
end
%% Try to save data to file
try
  fprintf('Saving bathymetry data to file: %s\n',fullfile(pwd,[project_name '_bathymetry.mat']))
  save(fullfile(pwd,[project_name '_bathymetry.mat']),'bathy','-v7.3'); % ,'-v7.3' is a compression flag
catch
  save(fullfile(pwd,[datestr(now,'yyyymmdd') '_bathymetry.mat']),'bathy','-v7.3'); % ,'-v7.3' is a compression flag
end
end %% FUNCTION get_bathymetry