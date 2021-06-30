function [plots] = UVP_select_paramets_to_plot(data,info)
%%  UVP_SELECT_PARAMETS_TO_PLOT
%
%  Syntax:
%    [plots] = UVP_select_paramets_to_plot(data,info)
%
%  Description:
%    Prompts the user to select which variables to plot. If user selects a
%    variable in NxM size, it asks again which size class or taxa to plot. 
%
%  Input:
%    data - ZOO, PAR, or CTD data structure
%    info - structure containing information about each field in "data"
%
%  Output:
%    plots - structure with fields of each variable to plot with
%    information on plot title and axis/color limits.
%
%  Dependencies:
%
%  Notes:
%    At the end of the script, it will print formatted text to the screen
%    so you don't have to run through this repeatedly. 
%    ... for example ...
%    *********************************************************************************
%    To hardcode these selections, copy and paste this into your UVP_workflow.m script
%    *********************************************************************************
%    plots = struct();
%    plots.abund.title = {'Phaeodendrida[#/m^3]'};
%    plots.abund.clims = [[0 0]];
%    *********************************************************************************
%
%  Author: Brita Irving <bkirving@alaska.edu>
%%
percent_limits = [0.01 99.5]; % percent limits used in plot color limits and axis limits. I.e. [0.5 99.5] = 0.5% to 99.5%
FontSize = 16;
%% Pull out list of available fields for plotting
fields = fieldnames(data);
idx_meta    = find(strcmp(fields,'SampledVolume'));
% Limit to data fields, do not allow user to select meta variables such as
% profile, volume, depth, etc.
data_fields = fields(idx_meta+1:end);
data_fields(contains(data_fields,{'rel_unc' 'datenum' 'header'})) = [];
%% Create Radio Button figure
try
  MP = get(0, 'MonitorPositions'); % Get monitor positions in pixel size
  if size(MP,1) > 1; nMP = 2; else; nMP = 1; end
  % Now change position to second monitor
  position = MP(nMP,:);
  position(1) = position(1)*1.1; % move left 10%
  position(2) = position(2)+60;  % move up 60 pixels
  position(3) = position(3)/2;   % shrink width
  position(4) = position(4)*0.9; % shrink heigh to 90%
  fig = uifigure('Position',position);
catch
  fig = uifigure;
end
fig.Name = 'Variable selection';
%% Create check boxes for all available data fields
lbl = uilabel(fig,'Text','Select which variables to plot');
lbl.Position(2) = fig.Position(4)-50;
lbl.Position(3) = fig.Position(3);
lbl.Position(4) = lbl.Position(4)*3;
lbl.FontSize    = 20; lbl.FontWeight  = 'bold';
% get height of figure so can space check boxes appropriately
ypos_start = fig.Position(4) - fig.Position(4)/10;
ypos_step  = floor(fig.Position(4)/numel(data_fields))/2;
% if only a few variables, change sizing so doesn't overtake screen
if numel(data_fields) < 5
  ypos_step = ypos_step/4;
end
% Loop through fields and create checkbox
for n = 1:numel(data_fields)
  checks(n) = uicheckbox(fig);
  try
    checks(n).Text = [pad(data_fields{n},20) ' | ' pad(info.(data_fields{n}).name,60)];
  catch
    if numel(info.(data_fields{n}).name) > 1
      checks(n).Text = [pad(data_fields{n},20) '    (You will be asked to select which fields next)'];
    else
      checks(n).Text = data_fields{n};
    end
  end
  checks(n).FontSize    = FontSize;
  checks(n).FontWeight  = 'bold';
  checks(n).Position(2) = ypos_start - ypos_step*n;
  checks(n).Position(3) = fig.Position(3)*0.9; % 90% of figure width
  checks(n).Position(4) = checks(n).Position(4)*2; % 90% of figure width
end

% Add checkbox for "FINISHED" so can move on
check_done = uicheckbox(fig);
check_done.Text = '* Done selecting variables * ';
check_done.FontSize    = 17;
check_done.FontWeight  = 'bold';
check_done.FontColor   = 'r';
check_done.Position(2) = ypos_start - ypos_step*(n+3);
check_done.Position(3) = fig.Position(3)*0.9; % 90% of figure width
check_done.Position(4) = checks(n).Position(4)*2; % 90% of figure width
% Wait until the "Done selectiong variables!" checkbox is checked
waitfor(check_done,'Value',1);

%% Pull out selected variable names
select_fields = data_fields([checks(1:end-1).Value]);

%% Create plotting characteristics based on selected variable
% "plots" struture will store the fields, their titles, and their color
% limits. 
plots = struct();
for nplot = 1:numel(select_fields)
  % Variable was checked
  field = select_fields{nplot};
  plots.(field) = struct();
  plots.(field).title = info.(field).name; % title for plotting
  
  % Change from char to cell so can pull out title by simple calling
  % plots.(field).title{nidx} rather than plots.(field).title(nidx,:).
  if ischar( info.(field).name)
    plots.(field).title = {plots.(field).title};
  end
  %% ONLY 1 Choice - e.g. Nx1 variable
  if numel(plots.(field).title) == 1
    plots.(field).clims = prctile(data.(field),percent_limits); % color limits: 0.5 to 99.5 %
    %% NxM variable
  else
    try % Create a new figure with specific position
      MP = get(0, 'MonitorPositions'); % Get monitor positions in pixel size
      if size(MP,1) > 1; nMP = 2; else; nMP = 1; end
      % Now change position to second monitor
      position = MP(nMP,:);
      position(1) = position(1)*1.1; % move left 10%
      position(2) = position(2)+60;  % move up 60 pixels
      position(3) = position(3)/2;   % shrink width
      position(4) = position(4)*0.9; % shrink heigh to 90%
      fig2 = uifigure('Position',position);
    catch % Create a new figure with default positioning
      fig2 = uifigure;
    end
    fig2.Name = [field 'variable selection'];
    % label the check boxes
    lbl = uilabel(fig2,'Text',['Select which ' field ' variables to plot']);
    lbl.Position(2) = fig2.Position(4)-50;lbl.Position(3) = fig2.Position(3);
    lbl.Position(4) = lbl.Position(4)*3;  lbl.FontSize    = 20; lbl.FontWeight  = 'bold';
    % mark half the fields so can split into two columns
    half = numel(info.(field).name)/2;
    % get height of figure so can space check boxes appropriately
    ypos_start = fig2.Position(4) - fig2.Position(4)/15;
    ypos_step  = floor(fig2.Position(4)*0.95/half);
    % Loop through fields and create checkbox
    for n = 1:numel(info.(field).name)
      if all(isnan(data.(field)(:,n)))
        continue
      end
      % Split into two columns
      
      if n < half
         xpos = fig2.Position(3)/10; 
         ypos = ypos_start - ypos_step*n;
      else
        xpos = fig2.Position(3)/2;              % start in the middle of the figure
        ypos = ypos_start - ypos_step*(n-half); % restart at top
      end
      checks2(n) = uicheckbox(fig2);
      checks2(n).Text = info.(field).name{n};
      checks2(n).FontSize    = 12; checks2(n).FontWeight  = 'bold';
      checks2(n).Position(1) = xpos; 
      checks2(n).Position(2) = ypos; 
      checks2(n).Position(3) = fig2.Position(3)*0.9;     % 90% of figure width
      checks2(n).Position(4) = checks2(n).Position(4)*2; %d ouble height
    end
    
    % Add option to select all
    check_all = uicheckbox(fig2);
    check_all.Text = ['Select ALL'];
    check_all.FontSize    = 16;check_all.FontWeight  = 'bold';check_all.FontColor   = 'k';
    check_all.Position(1) = lbl.Position(1); %fig2.Position(3)/2; %halfway along figure
    check_all.Position(2) = fig2.Position(4)*0.03; % 3% above the bottom
    check_all.Position(3) = fig2.Position(3)*0.9;  % 90% of figure width
    check_all.Position(4) = lbl.Position(4); %check_all.Position(4)*2; % double height
    
    % Add checkbox for "FINISHED" so can move on
    check_done2 = uicheckbox(fig2);
    check_done2.Text = ['* Done selecting ' field ' variables * '];
    check_done2.FontSize    = 17; check_done2.FontWeight  = 'bold'; check_done2.FontColor   = 'r';
    check_done2.Position(1) = lbl.Position(1);
    check_done2.Position(2) = 0.00;% very bottom
    check_done2.Position(3) = fig2.Position(3)*0.9; % 90% of figure width
    check_done2.Position(4) = lbl.Position(4);% check_done2.Position(4)*2; % double height
    % Wait until the "Done selectiong variables!" checkbox is checked
    waitfor(check_done2,'Value',1);
    
    % Catch if no variables selected
    if check_all.Value == 0 &&  all([checks2.Value] == 0)
      plots = rmfield(plots,field);
      continue % Skip to next variable
    end
    % Initialize indice array to catch nan cases
    rm_these_indices = [];
    % Initialize clims
    if check_all.Value == 0
      plots.(field).title = plots.(field).title([checks2.Value]);
      plots.(field).clims = nan(sum([checks2.Value]),2); % intialize
      data_idx = find([checks2.Value] == 1);
    else
      plots.(field).clims = nan(numel(plots.(field).title),2); % intialize
      data_idx = 1:numel(plots.(field).title);
    end
    % Loop through array and, if variable was checked, find plot limits
    for nvar = 1:numel(plots.(field).title)
      idx_here = data_idx(nvar);
      if all(isnan(data.(field)(:,idx_here)))
        % Catch if all data is NaN -- this could happen if did not delete
        % empty size bins in UVP_read_odv_ecotaxa_exported_par.m
        rm_these_indices = [rm_these_indices; nvar];
      else % contains real data
        plots.(field).clims(nvar,:) = prctile(data.(field)(:,idx_here),percent_limits);
      end
    end % Loop through data indices to pull out limits
    % Remove nan indicies, if any were found.
    plots.(field).title(rm_these_indices)   = [];
    plots.(field).clims(rm_these_indices,:) = [];
    delete(fig2); clearvars checks2 check_done2 check_all rm_these_indices
  end % IF Nx1 OR NxM VARIABLE/FIELD
end % LOOP THROUGH SELECTED VARIABLES/FIELDS
% Close uifigures
delete(fig);

%% Print selected variables to screen so can hardcode next time
fields = fieldnames(plots);
fprintf('*********************************************************************************\n')
fprintf('To hardcode these selections, copy and paste this into your UVP_workflow.m script\n')
fprintf('*********************************************************************************\n')
fprintf('plots = struct();\n')
for nf = 1:numel(fields)
  field = fields{nf};
  titles = [];
  clims  = [];
  for nt = 1:numel(plots.(field).title)
    titles = [titles '''' plots.(field).title{nt} ''';'];
    clims  = [clims '[' num2str(plots.(field).clims(nt,1)) ' ' num2str(plots.(field).clims(nt,2)) '];'];
  end
  if strcmp(titles(end),';'); titles(end) = []; end
  if strcmp(clims(end),';');  clims(end)  = []; end
  fprintf('plots.%s.title = {%s};\n',field,titles);
  fprintf('plots.%s.clims = [%s];\n',field,clims);
end
fprintf('*********************************************************************************\n')

end %% MAIN FUNCTION : UVP_SELECT_PARAMETS_TO_PLOT