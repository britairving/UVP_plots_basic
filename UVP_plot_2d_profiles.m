function UVP_plot_2d_profiles(data,info,plots)
%% FUNCTION UVP_PLOT_2D_PROFILES
%  ------------------------------------------------------------------------
%  Syntax: UVP_plot_2d_profiles(data,info,plots)
%
%  Description:
%    Plots data exported from ecotaxa AND UVP_process_options.m
%
%  Inputs:
%    opt  = structure with project options so don't need to edit code
%    directly. Input project information into UVP_[project]_options.m
%    uvp = data structure from reformat_UVP_process_to_ecotaxa_output.m
%           -----------------------OR-------------------------
%        = data structure from UVP_read_ecotaxa_exported.m
%
%  Dependencies: (in google drive & github repository)
%     utility folder - contains necessary dependent scripts
%
%  Notes:
%     Script uses subfunctions for organization.
%     Blue variables are global - accessible across all subfunctions
%     Plotting
%        set standard figure configurations in uvp_makefig.m
%        set standard print figure configurations in uvp_printfig.m
%        makefig_subplots(3,1) = subplot(1,3,1); subplot(1,3,2); subplot(1,3,3)
%        Some plotting scripts are subfunctions while others call external functions (i.e. plot_station_map.m and UVP_3d_plots.m)
%
%  Authors: Brita Irving <bkirving@alaska.edu>
%  ------------------------------------------------------------------------
%% Define all script options here
warning('off')

script_opt = struct();
script_opt.savefig             = 0;    % 0=no 1=yes -- saved to script_opt.savedir
script_opt.savedir             = opt.path.analysis;
script_opt.initial_config      = 0;    % stops after initial plot so can adjust parameters as you'd like, also pauses each transect iteration
script_opt.plot_in_TSspace     = 0;    % 1=plots variables defined in p.vars in TS diagram
script_opt.log_parabundance    = 0;    % 1=uses log colorbar for tot_par_abundance contourf plot
script_opt.yaxis_dens          = 0;    % 0=yaxis depth 1=DOES NOT LOOK GOOD, DON'T USE?? plot density as yaxis instead of depth
script_opt.ignore_profiles     = {};   % reads opt.ignore_profiles & skips selected profiles during transect plots
script_opt.show_profile_locs   = 1;    % 1=shows profile locations on contourf plots
script_opt.grid_by_time        = 1;    % 1=grids by time and overrides grid_by_lat
script_opt.grid_by_lat         = 0;    % 1=grids by latitude, 0=grids by longitude
script_opt.closefig            = 0;    % 0=no 1=yes
script_opt.zoo                 = 1;    % 1=zoo data read from ecotaxa file, 0=only par data
script_opt.botl                = 0;    % 1=botl data available
script_opt.ctd                 = 0;    % 1=ctd data available, 0=not available
script_opt.split_1000m         = 1;    % plots first 1000m as 50% of axis
script_opt.depth_bin           = 10;    % grid to x meters - for shallow cruises use small value
script_opt.limit_interpolation = 0;    % 1 = masks interpolation over vast distances, 0 = no mask applied - will interpolate over everything
script_opt.axcolor             = [0.9 0.9 0.9];  % axis color
script_opt.contourlevs         = 250;      % more levels = more detailed
script_opt.colormap            = 'lansey'; % lansey = lansey.m, 'cm' = lowval_colormap.m (highlights low values)
script_opt.colormap_profbysize = 'lines';  % different colormap because "particle profiles by size class" plot difficult to see
script_opt.loop_sections       = 1;    % plots individual sections - must define in opt.sections
script_opt.start_profile       = 1;    % Set to profile number to start plotting
% I.E. change script_opt.start_profile if on cruise
script_opt.plot_each_profile   = 0;    % plots individual profile statistics
script_opt.plot_transects_wmap = 0;    % plots multiple transects in single figure with map
script_opt.plot_variable_trans = 1;    % plots single variable transects
script_opt.plot_CSD_transects  = 0;    % plots all diameter bin CSDs (i.e. 27+ plots!)
script_opt.plot_CSD_comb_trans = 0;    % plots all diameter bin CSDs combined in 6 subplots
script_opt.plot_VSD_comb_trans = 0;    % plots all diameter bin volume combined in 6 subplots
script_opt.plot_3d_transects   = 0;    % plots 3D transects
%% Change x axis for section plots
% if strcmp(opt.cruise_name,'sn207_2018_s04p') % GRID BY LONGITUDE
%   script_opt.grid_by_time        = 0;    % 1=grids by time and overrides grid_by_lat
%   script_opt.grid_by_lat         = 0;    % 1=grids by latitude, 0=grids by longitude
% % elseif strcmp(opt.cruise_name,'sn207_2018_exports_np_sr1812')
% %   script_opt.grid_by_time        = 1;    % 1=grids by time and overrides grid_by_lat
% %   script_opt.grid_by_lat         = 0;    % 1=grids by latitude, 0=grids by longitude
% end

%% configure data structure if zoo data available
if isfield(uvp,'zoo')
  script_opt.zoo = 1;
  zoo = uvp.zoo;
  uvp = uvp.par;
else
  script_opt.zoo = 0;
end

%% Check data structure is acceptable format
if ~isfield(uvp,'concen')
  fprintf('Data not in correct format\n')
  fprintf('  If plotting data from UVP_Process_options.m then need to run it through "reformat_UVP_process_to_ecotaxa_output.m"\n')
end



%% Check transect/sections structure is acceptable format
% Definitions:
%   section  = range of profiles that define a specific part of the cruise
%   transect = old name for section but confusing because transect also
%              mean the contourf plots
if script_opt.loop_sections && isfield(opt,'transects') && ~isfield(opt,'sections')
  opt.sections = opt.transects;
  for is = 1:numel(opt.transects)
    pid1 = opt.transects{is};
    if is == numel(opt.transects)
      pid2 = meta.profileid{end};
    else
      pid2 = opt.transects{is+1};
    end
    % catch simply case where profileid's contain "_"
    if contains(pid1,'_') && contains(pid2,'_')
      rm_str = pid1(1:strfind(pid1,'_'));
      n1 = str2num(erase(pid1,rm_str));
      n2 = str2num(erase(pid2,rm_str)) - 1; % early definition just listed first profile. Now need list of all profiles in each section
    else
      rm_str = pid1(isletter(pid1));
      n1 = str2num(erase(pid1,rm_str));
      n2 = str2num(erase(pid2,rm_str)) - 1;
    end
    if isempty(n1) || isempty(n2)
      fprintf('** error here: need to change way sections are formatted--used for looping through and plotting sections/transects from cruise\n')
      keyboard
    end
    step1 = num2str(n1:n2,' %03d');
    step2 = strsplit(step1);
    step3 = strcat(rm_str, step2);
    opt.sections{is} = step3;
  end
end


%% Select which plots to visualize
% Plot Individual Profile data
% if script_opt.plot_each_profile
%   loop_profiles_plot_data;
% end

if script_opt.ctd
  plot_cruise_ctd
end
%
% Waterfall plot mips & maps
% plot_mips_maps_waterfall

% Waterfall plot total particle abundance
% plot_total_particle_abundance_waterfall

% Plot size distributions at each station
% plot_station_sizes

if script_opt.plot_in_TSspace
  plot_TSspace
end
% Plot Transect data
plot_transect_data;

% Plot 3D Transects
if script_opt.plot_3d_transects
  script_opt.p = p; % pass through variables
  UVP_3d_plots_sr_rachel(uvp, opt, meta, script_opt)
end

%% FUNCTION plot_station_sizes
  function plot_station_sizes
    if isnan(uvp.sizeinfo.szav(end))
      uvp.sizeinfo.szav(end) = 26;
    end
    minCSDn = min(min(uvp.CSDn));
    maxCSDn = max(max(uvp.CSDn));
    cmap_csd = cmapping(log(minCSDn:maxCSDn),script_opt.colormap,'continuous');
    ax = makefig_subplots(1,2);
    % now loop through profiles
    for ip = 1:script_opt.num_profiles
      prof = meta.profileid{ip};
      prng = find(contains(uvp.profile,prof));
      if isempty(prng)
        fprintf('Skipping profile %s\n',prof)
        continue
      end
      plot(ax(2),uvp.CSD(prng,:),uvp.Depth(prng),'x','MarkerSize',4,'Color','k')
      set(ax(2), 'XScale','log');
      title(ax(2),[strrep(opt.cruise_name,'_',' ') ' Profile: ' strrep(prof,'_','\_') ' CSD'])
      c = pcolor(ax(1),uvp.sizeinfo.szav,uvp.Depth(prng),uvp.CSDn(prng,:));
      colormap(ax(1),cmap_csd);
      clims = [minCSDn maxCSDn];
      c.FaceColor = 'interp'; c.EdgeColor = 'none'; cbar = colorbar(ax(1));
      cbar.Label.String = 'CSDn';
      ax(1).CLim = clims;
      cbar.Limits = clims;
      %contour(ax(1),uvp.sizeinfo.szav,uvp.Depth(prng),ncsd(prng,:),100)
      if script_opt.initial_config
        fprintf(' stopped here because set initial_config = 1, "dbcont" to continue\n')
        keyboard
      end
      if script_opt.savefig
        uvp_printfig(fullfile(script_opt.savedir,[opt.project_name '_' prof '_CSD']))
      end
    end
  end %% FUNCTION plot_station_sizes


%% FUNCTION plot_total_particle_abundance_waterfall
  function plot_total_particle_abundance_waterfall
    uvp_makefig;ax = gca;set(ax,'Color',script_opt.axcolor)
    title(ax,[strrep(opt.cruise_name,'_',' ') ' Total Particle Abundance'])
    xlabel(ax,'UVP Abundance [#/L]')
    ylabel(ax,'Depth [m]'); set(gca,'YDir','reverse');
    hold(ax,'on');
    int = 5;
    for ip = 1:script_opt.num_profiles
      prof = meta.profileid{ip};
      prng = find(contains(uvp.profile,prof));
      if isempty(prng)
        fprintf('Skipping profile %s\n',prof)
        continue
      end
      plot(ax(1),uvp.tot_par_abundance(prng)+(int*ip),uvp.Depth(prng),'-','Color','k')
    end
    grid(ax,'on')
    
    % highlight first 1000m
    if script_opt.split_1000m
      xlim1 = ax.XLim;
      % change limits and shrink positions
      ax.YLim = [1000 script_opt.maxdepth];
      ax.Position(4) = ax.Position(4)/2;
      ax.YAxis.Label.Units = 'normalized';
      ax.YAxis.Label.Position(2) = 1;
      % copy axes with data
      ax_top = copyobj(ax,gcf);
      % change limits and shrink positions
      ax_top.YLim = [0 1000]; ax_top.Position(2) = ax_top.Position(2) + ax_top.Position(4);
      ax_top.YLabel     = [];
      ax_top.XTickLabel = [];
      ax_top.XLabel     = [];
      % make sure xlimis consistent
      ax.XLim        = xlim1;
      ax_top.XLim       = xlim1;
    end
    % set axis labels to zero since not quantitative
    ax.XTickLabel = [];
    text(ax,0.9,0.1,['Total Particle Abundance offset by ' num2str(int) ' per profile'],'units','normalized','HorizontalAlignment','right')
    if script_opt.savefig
      uvp_printfig(fullfile(script_opt.savedir,[opt.project_name '_waterfall_total_par_abundance']))
    end
  end %% FUNCTION plot_total_particle_abundance_waterfall


%% FUNCTION plot_TS_space
  function plot_TSspace
    for iv = 1:numel(p.vars)
      uvp_makefig; ax= gca;
      fprintf(' Plotting %s in TS space\n',p.ptitle.(p.vars{iv}))
      plot_vars_TSspace(p.vars{iv},p.ptitle.(p.vars{iv}),p.clims.(p.vars{iv}));
      set(gca,'Color',[0.9 0.9 0.9])
      
      if script_opt.initial_config
        fprintf(' stopped here because set initial_config = 1, "dbcont" to continue\n')
        keyboard
      end
      if script_opt.savefig
        save_figname = fullfile(script_opt.savedir,[opt.project_name '_' p.stitle.(p.vars{iv}) '_TS']);
        if  script_opt.log_parabundance && contains(p.vars{iv},'tot_par_abundance')
          save_figname = [save_figname '_logcolorbar'];
        end
        uvp_printfig(save_figname)
      end
    end
  end

%% FUNCTION plot_cruise_ctd
  function plot_cruise_ctd
    num_axes = 5;
    ax = makefig_subplots(num_axes,1);
    
    % plot ctd data
    fnames = fieldnames(uvp);
    % salinity
    salt.fname = fnames{contains(fnames,{'practical_salinity','practicalSalinity'})};
    salt.name = 'practical salinity [psu]';
    % temperature
    temp.fname = fnames{contains(fnames,{'temperature'})};
    temp.name = 'temperature [degc]';
    % Oxygen
    idx = find(contains(fnames,{'oxygen','Oxygen'}));
    if ~isempty(idx)
      oxy.fname = fnames{idx};
      oxy.name = strrep(oxy.fname,'_',' ');
    end
    % pressure
    pres.fname = fnames{contains(fnames,{'pressure'})};
    pres.name = 'pressure [db]';
    % density
    idx = find(contains(fnames,{'density'}));
    if ~isempty(idx)
      dens.fname = fnames{idx};
      dens.data = uvp.(dens.fname);
      dens.name = strrep(dens.fname,'_',' ');
    else
      dens.data = sw_dens(uvp.(salt.fname),uvp.(temp.fname),uvp.(pres.fname));
      dens.name = 'density [kg/m^3]';
    end
    % chlorophyll
    idx = find(contains(fnames,{'chloro'}));
    if ~isempty(idx)
      chl.fname = fnames{idx};
      chl.data = uvp.(chl.fname);
      chl.name = strrep(chl.fname,'_',' ');
    end
    
    ylims = [0 script_opt.maxdepth];
    if script_opt.num_profiles > 40
      fprintf(' only plotting first 40 profiles, can change here if desired\n')
      maxnump = 40;
    else
      maxnump =  script_opt.num_profiles;
    end
    
    for ip = 1:maxnump
      prof = meta.profileid{ip};
      pidx = find(contains(uvp.profile,prof));
      % total particle abundance
      plot_against_pressure(ax(1),uvp.tot_par_abundance(pidx),uvp.Depth(pidx),ylims,clrs(ip,:))
      
      % temperature
      plot_against_pressure(ax(2),uvp.(temp.fname)(pidx),uvp.(pres.fname)(pidx),ylims,clrs(ip,:))
      
      % salinity
      plot_against_pressure(ax(3),uvp.(salt.fname)(pidx),uvp.(pres.fname)(pidx),ylims,clrs(ip,:))
      
      % density
      plot_against_pressure(ax(4),dens.data(pidx),uvp.(pres.fname)(pidx),ylims,clrs(ip,:))
      
      if exist('oxy','var') % oxygen
        plot_against_pressure(ax(5),uvp.(oxy.fname)(pidx),uvp.(pres.fname)(pidx),ylims,clrs(ip,:))
      elseif exist('chl','var')% chlorophyll
        plot_against_pressure(ax(5),uvp.(chl.fname)(pidx),uvp.(pres.fname)(pidx),ylims,clrs(ip,:))
      else
        delete(ax(5));
      end
    end
    ax(1).XLabel.String = 'UVP Abundance [#/L]';
    ax(2).XLabel.String = temp.name;
    ax(3).XLabel.String = salt.name;
    ax(4).XLabel.String = dens.name;
    if exist('oxy','var')
      ax(5).XLabel.String = oxy.name;
    elseif exist('chl','var')
      ax(5).XLabel.String = chl.name;
    end
    % remove ylabels from ctd data
    for iax = 2:num_axes
      ax(iax).YLabel.String = [];
      ax(iax).YTickLabel = [];
    end
    % title on middle plot
    if mod(num_axes,2) % odd number of axes
      title(ax((num_axes+1)./2),strrep(opt.cruise_name,'_','\_'))
    end
    if script_opt.savefig
      uvp_printfig(fullfile(script_opt.savedir,[opt.project_name '_ctdplot']))
    end
  end %% FUNCTION plot_cruise_ctd

%% FUNCTION loop_profiles_plot_data
  function loop_profiles_plot_data
    for it = script_opt.start_profile:script_opt.num_profiles
      ProfileID = meta.profileid{it};
      % Load BRU & DAT files
      fprintf(' Plotting UVP profile %s data\n',ProfileID)
      % pull out profile data
      prof = meta.profileid{it};
      prng = find(contains(uvp.profile,prof));
      if isempty(prng)
        fprintf('Skipping profile %s\n',prof)
        continue
      end
      % Colormap for different sizes
      if isfield(script_opt,'colormap_profbysize')
        if ischar(script_opt.colormap_profbysize)
          eval(['cmap_profbysize=flipud(colormap(' script_opt.colormap_profbysize '(numel(uvp.sizeinfo.name))));']); close
        else
          cmap_profbysize=flipud(colormap(script_opt.colormap_profbysize(numel(uvp.sizeinfo.name)))); close
        end
      else
        cmap_profbysize=colormap(cmap(numel(uvp.sizeinfo.name))); close
      end
      %% Plot profile statistics
      ylims = [0 script_opt.maxdepth];
      
      %% Plot number of images per depth bin, total particle concentration, mean grayscale, mean size
      if isfield(uvp,'grayscale')
        ax = makefig_subplots(3,1); set(gcf,'Name',opt.cruise_name);
        
        plot_against_pressure(ax(1),uvp.tot_par_abundance(prng),uvp.Depth(prng),ylims,'k');
        ax(1).XLabel.String = 'Total Particle Abundance [#/L]';
        
        plot_against_pressure(ax(2),uvp.grayscale(prng),uvp.Depth(prng),ylims,'k');
        ax(2).XLabel.String = 'Mean Greyscale';
        ax(2).YLabel = []; ax(2).YTickLabel = [];  ax(2).Title.String = [];
        
        plot_against_pressure(ax(3),uvp.meansize(prng),uvp.Depth(prng),ylims,'k');
        ax(3).XLabel.String = 'Mean Size [mm]';
        ax(3).YLabel = []; ax(3).YTickLabel = [];  ax(3).Title.String = [];
      else
        ax = makefig_subplots(1,1); set(gcf,'Name',opt.cruise_name);
        plot_against_pressure(ax(1),uvp.tot_par_abundance(prng),uvp.Depth(prng),ylims,'k');
        ax(1).XLabel.String = 'Total Particle Abundance [#/L]';
      end
      
      title(ax(1),['Profile ' strrep(ProfileID,'_','-')],'HorizontalAlignment','right')
      % highlight first 1000m
      if script_opt.split_1000m
        for ia = 1:numel(ax)
          a = ax(ia);
          xlim1 = a.XLim;
          % change limits and shrink positions
          a.YLim = [1000 script_opt.maxdepth]; a.Position(4) = a.Position(4)/2;
          a.YAxis.Label.Units = 'normalized';
          a.YAxis.Label.Position(2) = 1;
          % copy axes with data
          ax_top = copyobj(a,gcf);
          % change limits and shrink positions
          ax_top.YLim = [0 1000]; ax_top.Position(2) = ax_top.Position(2) + ax_top.Position(4);
          ax_top.YLabel     = [];
          ax_top.XTickLabel = [];
          ax_top.XLabel     = [];
          % make sure xlimis consistent
          a.XLim       = xlim1;
          ax_top.XLim  = xlim1;
        end
      end
      if script_opt.savefig
        uvp_printfig(fullfile(script_opt.savedir,[opt.project_name '_' ProfileID '_ParticleProfile']))
      end
      
      %% Plot Concentration Size Distribution (CSDn) at all depths within one station
      uvp_makefig;ax = gca; hold(ax,'on');
      i=0;
      dep  = uvp.Depth(prng);
      for k=1:50:numel(dep)
        i=i+1;
        CSDnk = uvp.CSDn(prng(k),:);
        szav  = uvp.sizeinfo.szav;
        iz    = find(CSDnk == 0);
        CSDnk(iz) = [];
        szav(iz)  = [];
        
        %fprintf('%s\n',[num2str(dep(k)) 'm'])
        % plot(szav,CSDnk,'Color',cmap_profbysize(i,:),'Marker','o','LineStyle','-','DisplayName',[num2str(round(dep(k))) 'm'])
        plot(szav,CSDnk,'Color','b','Marker','o','LineStyle','-','DisplayName',[num2str(round(dep(k))) 'm'])
      end
      legend(ax,'show');
      ax.YDir = 'normal';
      set(gca,'Color',script_opt.axcolor)
      set(gca, 'YScale', 'log', 'XScale','log','Xlim',[0.05,25],'YLim',[10^-3,10^4]);
      xlabel('Equivalent Spherical Diameter [mm]');
      ylabel('Concentration Size Distribution [#/L/mm]');
      hold off
      title({strrep(ProfileID,'_', ' ');' PSDs by depth'});
      set(gcf,'Name',[strrep(meta.cruise{1},'_',' ') ' PSDs by depth'])
      if script_opt.savefig
        uvp_printfig(fullfile(script_opt.savedir,[opt.project_name '_' ProfileID '_CSDn']))
      end
      
      %% Plot vertical profiles of different sizes at this station
      uvp_makefig;ax=gca;
      set(gcf,'Name',[strrep(meta.cruise{1},'_',' ') ' Particle profiles by size class'])
      hold(ax,'on');
      for kp = 1:script_opt.num_sizes
        if all(uvp.CSD(prng,kp) == 0)
          continue
        end
        %fprintf('%s\n',[num2str(uvp.sizeinfo.szav(kp)) 'mm'])
        plot(ax,smoothdata(uvp.CSD(prng,kp),5),uvp.Depth(prng),'x','Color',cmap_profbysize(kp,:),'MarkerSize',4,'LineWidth',2,'DisplayName',[num2str(uvp.sizeinfo.szav(kp)) 'mm'])
        %plot(ax,smoothdata(uvp.CSD(prng,kp),5),uvp.Depth(prng),'Color',cmap_profbysize(kp,:),'Marker','.','LineStyle','-','DisplayName',[num2str(uvp.sizeinfo.szav(kp)) 'mm'])
      end
      hl = legend(ax,'show');
      set(ax,'XScale','log','XLim',[10^-3,10^3])
      set(ax,'Color',script_opt.axcolor)
      ylabel(ax,'Depth [m]')
      title(ax,{[strrep(uvp.cruise{1},'_','\_') ' ' strrep(ProfileID,'_', ' ')];' Smoothed particle profiles by size class'});
      
      if script_opt.savefig
        uvp_printfig(fullfile(script_opt.savedir,[opt.project_name '_' ProfileID '_SizesProf']))
      end
      if script_opt.closefig
        close all
      end
      if script_opt.initial_config
        fprintf('  stopped here because set initial_config = 1, "dbcont" to continue\n')
        keyboard
        fprintf('  skipping similar plots for other profiles\n')
        break
      end
    end
  end %% FUNCTION loop_profiles_plot_data

%% FUNCTION plot_against_pressure
  function plot_against_pressure(ax,x,y,ylim,clr)
    if nargin < 5; clr = 'k'; end
    axes(ax); hold(ax,'on'); grid(ax,'on')% make sure proper axes is current
    plot(ax,x,y,'x','color',clr,'LineWidth',2,'MarkerSize',4);
    set(ax,'XAxisLocation','bottom','YDir','reverse');
    ylabel(ax,'Pressure [dbar]')
    if nargin > 3
      set(ax,'YLim',ylim)
    else
      set(ax,'YLim',[0 script_opt.maxdepth])
    end
    set(ax,'Color',script_opt.axcolor)
    %title(strrep(meta.cruise{1},'_',' '))
    title(ax, strrep(meta.cruise{1},'_',' '))
  end %% FUNCTION plot_against_pressure

%% FUNCTION plot_vars_TSspace    
  function plot_vars_TSspace(var,axtitle,clims) %plot_vars_TSspace(Xl,Yl,var_data,p.ptitle.(p.vars{iv}),p.clims.(p.vars{iv}));
    if script_opt.log_parabundance && contains(axtitle,'tot_par_abundance')
      cmap_numL = cmapping(log10(clims(1):clims(2)),'-lansey_rev','continuous');
      colormap(cmap_numL);
    else
      colormap(cmap);
    end
    
    S = uvp.practicalSalinity_psu_;
    T = uvp.temperature_degc_;
    P = uvp.pressure_db_;
    C = uvp.(var);
    axislim = [min(S) max(S) min(T) max(T)];
    
    ax = tsdiagram_scatter(S,T,P,C,axislim);
    axtitle = strrep(axtitle,'_','\_');
    title(ax,axtitle)
    
    hcbar = colorbar(ax);
    hcbar.Label.String = axtitle;
    if nargin > 2 && ~isempty(clims) % climits defined
      ax.CLim = clims;
      hcbar.Limits = clims;
    else % climits not passed to subfunction
      try
        % based on percentiles
        percentfiles = prctile(reshape(z,numel(z),1),[0.5 99.5]); % do not count 0.5% outliers
        clims(1) = floor(percentfiles(1)/0.25)*0.25;% floor(percentfiles(1));
        clims(2) = ceil(percentfiles(2)/0.25)*0.25;% ceil(percentfiles(2));
        ax.CLim = clims;
        hcbar.Limits = clims;
      catch
        % auto
      end
    end
  end %% FUNCTION plot_vars_TSspace

%% FUNCTION standard_contourf
  function ax = standard_contourf(x,y,z,axtitle,clims)
    colormap(cmap);
    ax = gca;
    if script_opt.log_parabundance && contains(axtitle,'tot_par_abundance')
      cmap_numL = cmapping(log10(clims(1):clims(2)),'-lansey_rev','continuous');
      colormap(cmap_numL);
    else
      colormap(cmap);
    end
    
    if nargin > 4 && ~isempty(clims) % climits defined
      clevels = linspace(clims(1),clims(2),script_opt.contourlevs);
      contourf(ax,x,y,z,clevels,'LineStyle','none','LineColor',[0 0 0],'Fill','on');
    else
      contourf(ax,x,y,z,script_opt.contourlevs,'LineStyle','none','LineColor',[0 0 0],'Fill','on');
    end
    if script_opt.yaxis_dens % sigmaT on yaxis
      ax.YLim = [script_opt.sig1 script_opt.sig2];
    end
    set(ax,'YDir','reverse','fontsize',14,'Layer','top','FontWeight','bold')
    if script_opt.grid_by_time
      xlabel(ax,'Date');
      datetick(ax,'x');
    elseif script_opt.grid_by_lat
      xlabel(ax,'Latitude');
    else
      xlabel(ax,'Longitude');
    end
    ylabel(ax,'Pressure (db)')
    axtitle = strrep(axtitle,'_','\_');
    title(ax,axtitle)
    hold(ax,'on'); grid(ax,'on');
    hcbar = colorbar(ax);
    hcbar.Label.String = axtitle;
    %ax.XLim = [min(meta.latitude) max(meta.latitude)];
    ax.XLim = [min(min(x)) max(max(x))];
    if nargin > 4 && ~isempty(clims) % climits defined
      ax.CLim = clims;
      hcbar.Limits = clims;
    else % climits not passed to subfunction
      try
        % based on percentiles
        percentfiles = prctile(reshape(z,numel(z),1),[0.5 99.5]); % do not count 0.5% outliers
        clims(1) = floor(percentfiles(1)/0.25)*0.25;% floor(percentfiles(1));
        clims(2) = ceil(percentfiles(2)/0.25)*0.25;% ceil(percentfiles(2));
        ax.CLim = clims;
        hcbar.Limits = clims;
      catch
        % auto
      end
    end
    % Mask out bathymetry
    if script_opt.grid_by_time
      deg = datenum(num2str(meta.filename),'yyyymmddHHMMSS');
    elseif script_opt.grid_by_lat
      try
        deg = convert_latlon_zooprocess_to_decimaldegrees(meta.latitude);
      catch
        fprintf('Need to add path where convert_latlon_zooprocess_to_decimaldegrees.m script is located\n')
        fprintf('This is found in [MPDL]\MATLAB_scripts\utility\n')
        keyboard
      end
    else
      lon = wrapTo360(meta.longitude);
      deg = convert_latlon_zooprocess_to_decimaldegrees(lon);
    end
    bxaxis   = [deg(1);deg;deg(end)];
    
    bdepth = [script_opt.maxdepth;meta.bottomdepth;script_opt.maxdepth];
    bdepth = smoothdata(bdepth,3);
    %fill(blat,bdepth,0.75*ones(1,3))
    karea( bxaxis, bdepth, 'basevalue', script_opt.maxdepth, 'facecolor', [1 1 1]*0.5, 'edgecolor', [1 1 1]*0.1);
    %     % Plot location of profile data
    %     [latg,depthg]=meshgrid(uvp.lat,uvp.depth_bins_mid);
    %     ddi=find(~isnan(reshape(uvp.C_total,[numel(uvp.C_total),1])));
    %     plot(ax,latg(ddi),depthg(ddi),'.k','MarkerSize',.02)
  end

%% FUNCTION format_position_labels
  function xt2 = format_position_labels(xt)
    xt2 = {};
    if script_opt.grid_by_lat
      for ix = 1:numel(xt)
        if str2double(xt{ix}) < 1 % check if south
          xt2{ix} = [num2str(-1*str2double(xt{ix})) '^{\circ}S'];
        else
          xt2{ix} = [xt{ix} '^{\circ}N'];
        end
        if str2double(xt{ix}) == 0
          xt2{ix} = 'Equator';
        end
      end
    else
      for ix = 1:numel(xt)
        if str2double(xt{ix}) < 1 % check if south
          xt2{ix} = [num2str(-1*str2double(xt{ix})) '^{\circ}E'];
        else
          xt2{ix} = [xt{ix} '^{\circ}E'];
        end
      end
    end
  end %% FUNCTION format_position_labels

%% FUNCTION zoom_to_sections
  function zoom_to_sections(zoom)
    top_axes = zoom.ax(end);
    % plot map so can see where sections are
    if ~isfield(zoom,'ax_map')
      zoom.ax.Position(1) = 0.18; % move axis to the right
      if script_opt.split_1000m
        zoom.ax_top.Position(1) = zoom.ax.Position(1);
        % shift colorbar too
        zoom.ax.Colorbar.Position(1) = zoom.ax.Position(1) + zoom.ax.Position(3) + 0.01;
      end
      zoom.ax_map = axes('Position',[0.015 0.75 0.12 0.1]); % map axes
      plot_station_map_small(opt,meta,zoom.ax_map);
      colormap(cmap);
      zoom.ax_map.YLabel.String = []; zoom.ax_map.XLabel.String = []; zoom.ax_map.Title.String = [];
      zoom.ax_map.XTickLabelRotation = 45;
    end
    % Loop through sections
    vec = zoom.vec;
    orig_title = zoom.ax(1).Title.String;
    for it = 1:numel(opt.sections)
      fprintf('   Plotting transect %d/%d\n',it,numel(opt.sections))
      % pull out profileID range
      % keyboard
      pidx = contains(transpose(vec.site),opt.sections{it});
      xrange = [min(vec.geo) max(vec.geo(pidx))]; % containts calls ctdxxx and station names. Try to make it call ctdxxx for both.
      % xrange = [min(vec.geo(pidx)) max(vec.geo(pidx))]; % original
      % adjust xaxis limits
      if size(zoom.ax,2) > 1
        for iv = 1:numel(p.vars)
          %zoom.ax(iv).Position(3) = 0.3; % not sure why used to do this...
          zoom.ax(iv).XLim = xrange;
          if script_opt.grid_by_time
            datetick(zoom.ax(iv),'x','keeplimits')
          end
          if script_opt.split_1000m
            top_axes =  zoom.ax_top(end);
            zoom.ax_top(iv).XLim = zoom.ax(iv).XLim;
            zoom.ax_top(iv).XTick = zoom.ax(iv).XTick;
          end
          
          if iv > 1
            zoom.ax(iv).XAxis.TickLabels = [];
          end
        end
      else
        zoom.ax.XLim = xrange;
        if script_opt.grid_by_time
          datetick(zoom.ax,'x','keeplimits')
        end
        if script_opt.split_1000m
          zoom.ax_top.XLim = zoom.ax.XLim;
          zoom.ax_top.XTick = zoom.ax.XTick;
          zoom.ax_top.Title.String = zoom.ax.Title.String;
          zoom.ax.Title.Visible = 'off';
        end
        
      end
      top_axes.Title.String = { orig_title; ...
        ['Section ',num2str(it),'/',num2str(numel(opt.sections)) ': Profiles ' strjoin(opt.sections{it})]};
      
      % highlight section on map, if there is a map in the figure
      if isfield(zoom,'ax_map')
        axes(zoom.ax_map);
        m_plot(vec.lon(pidx),vec.lat(pidx),'og','MarkerFaceColor','r','MarkerSize',8,'clipping','on');
        title(zoom.ax_map,{strrep(uvp.cruise{1},'_','\_'); ['Section ',num2str(it),'/',num2str(numel(opt.sections))]})
      end
      % Save the figure
      if isfield(zoom,'save_figname')
        uvp_printfig([zoom.save_figname '_transect' num2str(it)])
      end
      
      % pauses here so can look at figures and get parameters set as desired
      % before looping through each variable/transect
      if script_opt.initial_config
        fprintf('  stopped here because initial_config == 1 ... "dbcont" to continue\n')
        keyboard
      end
      % If map exsits - erase section highlight
      if isfield(zoom,'ax_map')
        m_plot(vec.lon(pidx),vec.lat(pidx),'og','MarkerFaceColor','k','MarkerSize',8,'clipping','on');
        title(zoom.ax_map,'')
      end
    end
  end %% FUNCTION zoom_to_sections
end