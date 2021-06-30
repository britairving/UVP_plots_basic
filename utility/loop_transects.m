    % zoom in on transects
    if script_opt.loop_transects && isfield(opt,'transects') && ~isempty(opt.transects)
      set(ax,'YLim',[0 script_opt.maxdepth]);
      for it = 1:numel(opt.transects)
        fprintf('   Plotting transect %d/%d\n',it,numel(opt.transects))
        tidx1 = find(contains(meta.profileid,opt.transects{it}));
        if it == numel(opt.transects)
          tidx2 = size(meta.profileid,1);
        else
          tidx2 = find(contains(meta.profileid,opt.transects{it + 1})) - 1;
        end
        [tlat1,tlon1] = convert_latlon_zooprocess_to_decimaldegrees(meta.latitude(tidx1), meta.longitude(tidx1));
        [tlat2,tlon2] = convert_latlon_zooprocess_to_decimaldegrees(meta.latitude(tidx2), meta.longitude(tidx2));
        vidx1 = find(vec.lat == tlat1 & vec.lon == tlon1,1,'first');
        vidx2 = find(vec.lat == tlat2 & vec.lon == tlon2,1,'last');
        ax.XLim = [vec.lat(vidx1) vec.lat(vidx2)];
        if script_opt.savefig
          uvp_printfig([save_figname '_transect' num2str(it)])
        end
        % pauses here so can look at figures and get parameters set as desired
        % before looping through each variable/transect
        if script_opt.initial_config
          fprintf('  stopped here because initial_config == 1 ... "dbcont" to continue\n')
          keyboard
          fprintf(' skipping similar plots for other variables\n')
          break
        end
      end
      % Close figure
      if script_opt.closefig
        close(gcf)
      end
      clearvars var_data Xlat Ydepth
    end %% LOOP THROUGH TRANSECTS