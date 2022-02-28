function [mips,maps] = UVP_select_abundance_size_limits(par_info)
%%  UVP_SELECT_ABUNDANCE_SIZE_LIMITS
%
%  Syntax:
%    [plots] = UVP_select_abundance_size_limits(par_info)
%
%  Description:
%    Prompts the user to select which size bins to use for MIPs and MAPs
%    calculations. Although the general definition is easy, the size ranges
%    are somewhat arbitrary.
%    MIPS = microscopic particles, or small micrometric particles
%    MAPS = macroscopic particles, or large macroscopic particles
%
%  Input:
%    par_info - structure containing information on PAR data (read in with
%    UVP_read_odv_ecotaxa_exported_par.m)
%
%  Output:
%    idx - structure containing indicies to calculate mips and maps.
%
%  Dependencies:
%
%  Notes:
%
%  Author: Brita Irving <bkirving@alaska.edu>
%%
FontSize = 16;
%% Pull out list of available fields for plotting
size_bins = par_info.size_bins.strnam;

%% Create Radio Button figure for MIPs size bin selections
mips_done = 0;
while ~mips_done
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
  lbl = uilabel(fig,'Text',{'Select size bins to use for MIPs calculation'; '   MIPs = microscopic particles, or small micrometric particles'});
  lbl.Position(2) = fig.Position(4)-50;
  lbl.Position(3) = fig.Position(3);
  lbl.Position(4) = lbl.Position(4)*3;
  lbl.FontSize    = 20; lbl.FontWeight  = 'bold';
  % get height of figure so can space check boxes appropriately
  ypos_start = fig.Position(4) - fig.Position(4)/10;
  ypos_step  = floor(fig.Position(4)/numel(size_bins))/2;
  % if only a few variables, change sizing so doesn't overtake screen
  if numel(size_bins) < 5
    ypos_step = ypos_step/4;
  end
  % Loop through fields and create checkbox
  for n = 1:numel(size_bins)
    checks(n) = uicheckbox(fig);
    checks(n).Text = pad(size_bins{n},60);
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
  
  %% Store MIPs definitions
  mips.idx_size_bin = find([checks.Value]);
  if mips.idx_size_bin == 2
    % if only 1 size bin selected, assume that's correct and just include
    % that single size bin... i.e. 161-203um
    mips.idx_size_bin = [mips.idx_size_bin; mips.idx_size_bin];
  end
  if numel(mips.idx_size_bin) ~= 2
    fprintf('Must select start AND end bin sizes.. try again\n')
    mips_done = 0;
    close(fig);
    continue
  end
  text_pos = check_done.Position;
  text_pos(2) = text_pos(2) - 100;
  txt_check = uilabel(fig,'Text','**Verify choice in command window**','Position',text_pos,'FontSize',16,'FontColor','r','FontWeight','bold');
  mips.name = strjoin(size_bins( [checks.Value]),' to ');
  fprintf(' MIPs abundance will range from: %s\n', mips.name)
  mips_chc = input(' Is that correct? <1/0> ');
  switch mips_chc
    case 0
      mips_done = 0;
    case 1
      mips_done = 1;
    otherwise % Default to yes
      mips_done = 1;
  end
  close(fig);
end %% WHILE ~MIPS_DONE
%% Create Radio Button figure for MAPs size bin selections
maps_done = 0;
while ~maps_done
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
  lbl = uilabel(fig,'Text',{['Select size bins to use for MAPs calculation  (Note MIPs = ' mips.name ')']; '   MAPs = macroscopic particles, or large macroscopic particles'});
  lbl.Position(2) = fig.Position(4)-50;
  lbl.Position(3) = fig.Position(3);
  lbl.Position(4) = lbl.Position(4)*3;
  lbl.FontSize    = 20; lbl.FontWeight  = 'bold';
  % get height of figure so can space check boxes appropriately
  ypos_start = fig.Position(4) - fig.Position(4)/10;
  ypos_step  = floor(fig.Position(4)/numel(size_bins))/2;
  % if only a few variables, change sizing so doesn't overtake screen
  if numel(size_bins) < 5
    ypos_step = ypos_step/4;
  end
  % Loop through fields and create checkbox
  for n = 1:numel(size_bins)
    checks(n) = uicheckbox(fig);
    checks(n).Text = pad(size_bins{n},60);
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
  
  %% Store MIPs definitions
  maps.idx_size_bin = find([checks.Value]);
  if numel(maps.idx_size_bin) ~= 2
    fprintf('Must select start AND end bin sizes.. try again\n')
    maps_done = 0;
    close(fig);
    continue
  end
  maps.name = strjoin(size_bins( [checks.Value]),' to ');
  text_pos = check_done.Position;
  text_pos(2) = text_pos(2) - 100;
  txt_check = uilabel(fig,'Text','**Verify choice in command window**','Position',text_pos,'FontSize',16,'FontColor','r','FontWeight','bold');
  fprintf(' MAPs abundance will range from: %s\n', maps.name)
  maps_chc = input(' Is that correct? <1/0> ');
  switch maps_chc
    case 0
      maps_done = 0;
    case 1
      maps_done = 1;
    otherwise % Default to yes
      maps_done = 1;
  end
  close(fig);
end
end %% MAIN FUNCTION : UVP_SELECT_PARAMETS_TO_PLOT