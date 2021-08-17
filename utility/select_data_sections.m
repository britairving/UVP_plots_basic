function idx_section = select_data_sections(opt)
%% function select_data_sections
%
%  Syntax:
%    select_data_sections(opt)
%
%  Description:
%    Walks through loading sections
%
%  Inputs:
%    opt | structure containing 'lat' 'lon' and 'time'
%
%  Author: Brita Irving <bkirving@alaska.edu>
%% Choose how to split up sections
fprintf('\nChoose how to split up sections\n')
fprintf(' <1> Number of days (default)\n')
fprintf(' <2> Click points on a map\n')
section_chc2 = input('Enter section choice: ');
%% OPTION 1 | Number of days
if isempty(section_chc2) || section_chc2 == 1
  fprintf('\nTotal number of days in this cruise: %d\n',round(max(opt.time)-min(opt.time)));
  done = 0;
  while ~done
    day_choice = input('Enter number of days for each section: ');
    if isnumeric(day_choice) && day_choice > 0
      done = 1;
    else
      fprintf('the number of days for each section must be a number\n')
    end
  end % WHILE done selection number of days
  % Initialize
  nsection = 1; % section number
  tstart = min(opt.time); % start time if sectuib
  idx_section = [];       % matrix to store section indicies
  done_day_intervals = 0; % done flag
  while ~done_day_intervals
    rng = find(opt.time >= tstart & opt.time < tstart+day_choice);
    if isempty(rng) && tstart > opt.time(end)
      done_day_intervals = 1;
    elseif  isempty(rng) && tstart < opt.time(end)
      % Update start date for next iteration
      tstart = tstart+day_choice;
    else
      idx_section(nsection,1:2) = [rng(1) rng(end)];
      % Update start date for next iteration
      tstart = tstart+day_choice;
      nsection   = nsection + 1;
    end
  end
  
  %% OPTION 2 | Click points on a map
else
  
  makefig;ax1 = subplot(4,1,1:2); ax2= subplot(4,1,3:4);  %ax3= subplot(4,1,4);
  scatter(ax1,opt.time,opt.lat,50,opt.time,'filled');
  datetick(ax1,'x','mm/dd'); grid(ax1,'on'); ylabel(ax1,'latitude');
  cb = colorbar(ax1); cb.TickLabels = datestr(cb.Ticks,'mm/dd');
  
  % Plot time vs longitude
  scatter(ax2,opt.time,opt.lon,50,opt.time,'filled'); ax2.Position(3) = ax1.Position(3);
  datetick(ax2,'x','mm/dd'); ylabel(ax2,'longitude');grid(ax2,'on');
  
  ht1 = text(ax1,0.01,0.97,'Select section boundaries by clicking on graph','units','normalized','Color','r','FontWeight','bold');
  ht2 = text(ax1,0.01,0.90,'Hit RETURN <enter> when you are finished','units','normalized','Color','r','FontWeight','bold');
  
  fprintf('  Click on graph where you want section boundaries\n');
  fprintf('  Hit ENTER <return> when you are finished\n');
  % allows you to select an unlimited number of points until you press the Return key.
  [x,~] = ginput;
  
  idx_section = nan(numel(x),2);
  rm_this = [];
  for nx = 1:numel(x)
    
    if nx == numel(x)
      rng1 = find(opt.time >= x(nx),1);
      if isempty(rng1)
        rm_this = [rm_this; nx];
        continue
      else
        rng = [rng1 numel(opt.time)];
      end
    else
      rng = find(opt.time >= x(nx) & opt.time < x(nx+1));
    end
    
    % save section indices
    if ~isempty(rng)
      idx_section(nx,:) = [rng(1) rng(end)];
    else
      rm_this = [rm_this; nx];
    end
    % Remove empty sections
    %bad = find(isnan(idx_section(:,1)) | isnan(idx_section(:,2)));
    %rm_this = [rm_this; bad];
    %idx_section(rm_this,:) = [];
  end % Loop through number of "clicks" or sections
  bad = find(isnan(idx_section(:,1)) | isnan(idx_section(:,2)));
  if ~isempty(bad)
    idx_section(bad,:) = [];
  end
  
end % CHOICE ON HOW TO SELECT SECTIONS



end %% MAIN SCRIPT
