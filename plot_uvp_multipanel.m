function save_figname = plot_uvp_multipanel(par,par_info,plot_vars,options)
%% FUNCTION PLOT_UVP_MULTIPANEL
%  Syntax:
%    SAVE_FIGNAME = PLOT_UVP_MULTIPANEL(PAR,PAR_INFO)
%  
%  Description:
%    Creates a plot with 3 axis where Depth is on the y-axis and shows
%    total particle abundance, total particle biovolume, and slope of PSD
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Check all necessary fields are available
if nargin < 3 % Set Default variables
  plot_vars = {'tot_par_abundance' 'tot_par_biovolume' 'slope_b'};
end
if ~all( isfield(par,['Depth' plot_vars]) )
  error('"par" does not have the necessary fields!')
end

%% Defaults
cfg.split_watercolumn = 1;
prcntile_depths = prctile(par.Depth,[1 95]);      % Get percentiles in case there is a single really deep profile
cfg.maxdepth    = ceil(prcntile_depths(2)/50)*50; % Round maximum depth to next 50m
%% Initialize figure
ax = makefig_subplots(numel(plot_vars),1);
% ax = flipud(ax);
pause(0.2); % pause to matlab a chance to catch up.
fig = gcf;  fig.Name = 'UVP profiles';
% Get list of unique profiles
profiles = unique(par.profile);
% Create a table so can store each profile's data and show the mean profile
data = table();
data.depth = unique(par.Depth);
%% Loop through fields to plot
for nax = 1:numel(plot_vars)
  a = ax(nax);
  var = plot_vars{nax};
  % Title
  switch var
    case {'tot_par_abundance' 'tot_par_abundance2'}
      title_str = 'Total Particle Abundance';
    case {'tot_par_biovolume' 'tot_par_biovolume2'}
      title_str = 'Total Particle Volume';
    case 'slope_b'
      title_str = 'Slope of PSD';
    case 'meansize'
      title_str = 'Particle mean size [mm]';
    otherwise
      title_str = strrep(par_info.(var).name,'_','\_');
  end
  % Show the project, if possible
  if isfield(par,'header') && nax == 2
    title_str = {strrep(par.header,'_','\_'); title_str};
  end
  title(a,title_str);
  a.XLabel.String =  strrep(par_info.(var).name,'_','\_');
  
  % Initialize table so can plot the average profle 
  data.(var) = nan(numel(data.depth),numel(profiles));
  
  % Loop through profiles
  for np = 1:numel(profiles)
    idx_profile = strcmp(par.profile,profiles{np});
    smo_dat = smooth(par.(var)(idx_profile),3,'moving');
    
    plot_against_depth(a,smo_dat,par.Depth(idx_profile),[0 cfg.maxdepth]);
    %a.XAxisLocation = 'top';
    % Pull out profile data at available depths so can calculate average
    data.(var)(ismember(data.depth,par.Depth(idx_profile)),np) = smo_dat;
  end
  
  % Plot average profile
  havg = plot_against_depth(a,nanmean(data.(var),2),data.depth,[0 cfg.maxdepth],'g');
  havg.LineWidth = 6;
  % Get percentiles of data so can adjust x axis accordingly
  percentiles = prctile(par.(var),[0.001 99.5]);
  a.XLim(2) = percentiles(2); %max(nanmean(data.(var),2));
  
  % remove depth labels
  if nax == numel(plot_vars)
    a.YAxisLocation = 'right';
  elseif nax > 1
    a.YTickLabel    = [];
    a.YLabel.String = [];
  end
  % highlight top of the watercolumn
  if cfg.split_watercolumn && cfg.maxdepth < 1000
    % Pass output ax_top back so can change characteristics as necessary
    ax_top = plot_split_axes(a);
  end
  % Add annotation for slope of PSD 
  if strcmp(var,'slope_b')
    apos = get(a,'Position');
    an = annotation('doublearrow',[apos(1) apos(1)+apos(3)],[0.3 0.3],'Color','w');
    an.LineWidth = 4.0;
    % Put text somewhere near the middle bottom of axes
    text(a,0.05,0.25,{'Steep slope';  '"smaller"'},'Units','normalized','Color','w','FontWeight','bold');
    text(a,0.7, 0.25,{'Shallow slope';'"Larger"'}, 'Units','normalized','Color','w','FontWeight','bold');
  end
end %% FIELDS TO PLOT
%% Save figure
if options.savefig
  % Generate filename to save figure to
  save_figname = [options.project '_overview_of_particle_profiles'];
  standard_printfig(save_figname);
end

end %% MAIN FUNCTION
