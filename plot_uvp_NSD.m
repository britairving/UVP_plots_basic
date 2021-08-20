function plot_uvp_NSD(par,par_info,options)
%% FUNCTION plot_uvp_NSD
%  Syntax:
%    PLOT_UVP_NSD(DATA,DATA_INFO,SAVE_CHOICE)
%
%  Description:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Debug mode
% close all
dbstop if error
xscale = 'log'; % 'log' or 'linear'
%% Initialize figure
makefig;
pause(2); % pause to matlab a chance to catch up.
fig = gcf;

ax = gca;

%% Pull out unique profiles and number of size bins
profiles  = unique(par.profile);
size_clrs = jet(size(par.NSD,2));
variables = {'NSD' 'DNSD' 'VSD' 'DVSD'};
var_title = {'Number size distribution' 'Differential number size distribution' ...
             'Volume size distribution' 'Differential volume size distribution'};
for nv = 1:numel(variables)
  var = variables{nv};
  % Skip if field does not exist
  if ~isfield(par,var)
    continue
  end
  fig.Name = var_title{nv};
  %% Loop through profiles
  for np = 1:numel(profiles)
    try
      idx_profile = strcmp(par.profile,profiles{np});
      title(ax,['Profile: ' strrep(profiles{np},'_','\_') ' | ' var_title{nv}])
      hold(ax,'on');
      ax.YLabel.String = 'Depth [m]';
    catch
      keyboard
    end
    %% Loop through size bins
    for nsz = 1:size(par.(var),2) % Loop
      % Plot profile data
      plot(ax,par.(var)(idx_profile,nsz),par.Depth(idx_profile),'Color',size_clrs(nsz,:),...
        'Marker','o','MarkerSize',3,...
        'LineWidth',2,'LineStyle','-','DisplayName',par_info.(var).size_bin_st{nsz})
    end %% SIZE BINS
    % set xaxis to log or linear scale
    ax.XScale = xscale;
    grid(ax,'on');
    ax.XLabel.String = [var_title{nv} ' [' par_info.(var).unit{1} ']'];
    % show legend
    try
      hl = legend(ax,'show');
      hl.FontSize = 14;
      hl.Location = 'eastoutside';
    catch
      axes(ax);
      legend;
      ax.Legend.FontSize = 14;
      ax.Legend.Location = 'eastoutside';
      % not working on Matlab MAC
    end
    %% Save figure
    if options.savefig
      standard_printfig([options.project '_' var '_' profiles{np}]);
    end
    pause(0.05); % wait half a second for matlab to catch up
    cla(ax);
  end %% PROFILES
  
end %% DATA TYPES (NSD, VSD, & differential forms)

end %% MAIN FUNCTION
