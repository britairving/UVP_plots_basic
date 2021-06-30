function ax = makefig_subplots(subplotsx,subplotsy)
% make standard figure
% taken from http://p-martineau.com/perfect-subplot-in-matlab/
n = 1;
while ishandle(n)
  n = n + 1;
end
fig = figure(n);
pause(0.1);
% First, set units to normalized so can set screen size
fig.Units = 'normalized';
% change size to just a bit smaller than screen size
fig.OuterPosition = [0.02 0.04 0.90 0.90];
% Next, move figure to second monitor if possible
% If the get(0, 'MonitorPositions') function worked, and there are 2
% monitors, position the figure so it is on the second screen.
try
  % Get monitor positions in pixel size
  MP = get(0, 'MonitorPositions');
  if size(MP,1) > 1
    % Change units to Pixels because that is what get(0,[]) function returns
    fig.Units = 'pixels';
    % Now change position to second monitor
    fig.Position(1) = MP(2,1)+5;
    pause(0.2); % for some reason.. need to pause to change screen
    if ~isequal(fig.Position(1),MP(2,1)+5)
      fig.Position(1) = MP(2,1)+5;
    end
  end
catch
  % plot will stay in first screen
end

% Default # of x and y axes
if nargin == 0
  subplotsx = 3;
  subplotsy = 5;
end
leftedge=1.2;
rightedge=1.0;
topedge=1.5;
bottomedge=1.8;

spacex=0.4;
spacey=0.3;%1.0;%0.3;

fig.Units  = 'Inches';
plotwidth  = fig.Position(3);
plotheight = fig.Position(4)*2;

% Generate axes position
sub_pos = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
%loop to create axes
for i=1:subplotsx
  for ii=1:subplotsy
    ax(i,ii)=axes('position',sub_pos{i,ii},'Box','on','Layer','top');%,'XGrid','off','XMinorGrid','off','FontSize',fontsize,);
  end
end

fig.PaperPositionMode = 'auto';
fig.PaperOrientation  = 'portrait';
set(fig, 'color', 'w')
set(fig,'DefaultAxesLineStyleOrder','-|--|-.|:');
tomLnCmap = [ 0 0 1; 1 0 0; 1 0.6 0; 0.85 0 0.95; 0 .95 0.95; 0.8 0.8 0; 0.8 0.5 0; 0.5 0 0; 0 0 0.5; 0 0.65 0.65];
%               b;     r;   orange;      ~m;        ~bc;          ~y;    dorange;      dr;     db;       dc;
% no greens as Tom is R-G color blind.
set(fig,'DefaultAxesColorOrder',tomLnCmap);
% axes defaults
set(fig,'DefaultTextFontname',    'Times');
set(fig,'defaultTextFontSize',    16);
set(fig,'DefaultTextFontName',    'Times');
set(fig,'DefaultAxesFontName',    'Times');
set(fig,'DefaultAxesFontSize',    16);
set(fig,'DefaultAxesFontWeight',  'Bold');
set(fig,'DefaultAxesColor',[0.75 0.75 0.75])
end