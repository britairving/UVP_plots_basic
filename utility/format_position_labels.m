function xt2 = format_position_labels(xt,geo_choice)
%% FUNCTION format_position_labels
%
%  Author: Brita K Irving  <bkirving@alaska.edu>
%%
xt2 = {};
switch geo_choice
  %% LATITUDE
  case {'lat' 'latitude' 'Latitude'}
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
  %% LONGITUDE
  case {'lon' 'longitude' 'Longitude'}
    for ix = 1:numel(xt)
      if str2double(xt{ix}) < 1 % check if west
        xt2{ix} = [num2str(-1*str2double(xt{ix})) '^{\circ}W'];
      else
        xt2{ix} = [xt{ix} '^{\circ}E'];
      end
    end
end %% LATITUDE OR LONGITUDE
end %% FUNCTION format_position_labels

