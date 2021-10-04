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
fprintf(' <3> Specific profiles by entering profile names\n')
fprintf(' <4> Read profiles from file\n')
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
elseif section_chc2 == 2
  % random symbols
  symb = {'<' 'h' 'o' 's' 'p' '>' 'd'};
  symb = repmat(symb,1,20); % repeat a bunch of times incase many sections
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
    % Visualize section chosen
    try
      hold(ax1,'on');hold(ax2,'on');
      % Plot time vs latitude
      plot(ax1,opt.time(idx_section(nx,1):idx_section(nx,2)),opt.lat(idx_section(nx,1):idx_section(nx,2)),symb{nx},'Color','k','LineWidth',1,'MarkerSize',15);
      % Plot time vs longitud
      plot(ax2,opt.time(idx_section(nx,1):idx_section(nx,2)),opt.lon(idx_section(nx,1):idx_section(nx,2)),symb{nx},'Color','k','LineWidth',1,'MarkerSize',15);
    end
  end % Loop through number of "clicks" or sections
  bad = find(isnan(idx_section(:,1)) | isnan(idx_section(:,2)));
  if ~isempty(bad)
    idx_section(bad,:) = [];
  end
%% OPTION 3 | Specific profiles
elseif section_chc2 == 3
  % Initialize section indice matrix
  idx_section = nan(100,2);
  
  [unique_profiles,~] = unique(opt.profile,'stable');
  fprintf('Available profiles\n')
  fprintf('%s\n',strjoin(unique_profiles,','))
  fprintf('When you are done selecting sections, use <return> to exit loop\n')
  done_selecting_profile_sections = 0;
  section_num = 1;
  
  while ~done_selecting_profile_sections
    fprintf('Section #%d profile range\n',section_num);
    chc_p1 = input('  profile1: ','s');
    chc_p2 = input('  profile2: ','s');
    if ismember(chc_p1,unique_profiles) && ismember(chc_p2,unique_profiles) 
      idx_section(section_num,1) = find(strcmp(opt.profile,chc_p1),1);
      idx_section(section_num,2) = find(strcmp(opt.profile,chc_p2),1,'last');
      section_num = section_num + 1;
    elseif isempty(chc_p1) && isempty(chc_p2)
      done_selecting_profile_sections = 1;
      rm_placeholders = isnan(idx_section(:,1)) & isnan(idx_section(:,2));
      idx_section(rm_placeholders,:) = [];
    else
      fprintf('profiles not found, try again\n')
    end
    
  end
  %% OPTION 4 | Read profiles from file
elseif section_chc2 == 4
  if isfield(opt,'filename_sectionsbyprofiles') && exist(opt.filename_sectionsbyprofiles,'file')
    section_file = opt.filename_sectionsbyprofiles;
  else
    fprintf('** File must be formatted as follows ** \n')
    fprintf('"section_p1","section_p2"\n')
    fprintf('section1_p1,section1_p2\n')
    fprintf('section2_p1,section2_p2\n')
    fprintf('section3_p1,section3_p2\n')
    fprintf('section4_p1,section4_p2\n')
    fprintf('...\n')
    section_file = input('Enter filename: ','s');
    if ~exist(section_file,'file')
      fprintf('tried to find file "%s", but unsuccessful\n',section_file)
      error('File does not exist...')
    end
  end
  sections = readtable(section_file,'FileType','text','HeaderLines',0);
  if ~ismember('section_p1',sections.Properties.VariableNames) || ~ismember('section_p2',sections.Properties.VariableNames)
    fprintf('file "%s" must have column names "section_p1" and "section_p2"\n',section_file)
    error('File does not have proper column names...')
  end
  num_sections = size(sections,1);
  idx_section = nan(num_sections,2);
  for nsect = 1:num_sections
    idx_section(nsect,1) = find(strcmp(opt.profile,sections.section_p1(nsect)),1);
    idx_section(nsect,2) = find(strcmp(opt.profile,sections.section_p2(nsect)),1);
  end
end % CHOICE ON HOW TO SELECT SECTIONS

end %% MAIN SCRIPT
