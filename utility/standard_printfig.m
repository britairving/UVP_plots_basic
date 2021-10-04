function standard_printfig(figname)
%% function standard_printfig(figname)
% print figure in standard way
% author: Brita Irving <bkirving@alaska.edu>
%% Automatically format figure name
figname = strrep(figname,' ','_');
figname = erase(figname,'[#/L]');
figname = erase(figname,'[mm^3/L]');
figname = erase(figname,'[#/m^3]');
figname = strrep(figname,'(','_');
figname = strrep(figname,')','_');
figname = strrep(figname,'\','per');
figname = strrep(figname,'µ','u');
figname = strrep(figname,'__','_');
if strcmp(figname(end),'_')
  figname(end) = [];
end
%% Save figure
fprintf(' saving figure to %s\n',figname)
image_type       = '-djpeg';
image_resolution = '-r250';
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(fig, 'InvertHardcopy', 'off')
print(figname,image_type,image_resolution)
end %% FUNCTION STANDARD_PRINTFIG(FIGNAME)