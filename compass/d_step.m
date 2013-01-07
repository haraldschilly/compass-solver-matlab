%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d_step: A step kind in the NMS procedure in COMPASS.  
% Copyright (C) 2012 Stefan Schmelzer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This file is part of COMPASS.

% COMPASS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% COMPASS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
  
% You should have received a copy of the GNU General Public License
% along with COMPASS.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%
% Comments on this file: 

% A d_step is taken if the distance criterion is satisfied, if d_count
% is not over its limit, and if the merit function value is not too high
% above Ref.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


z_k = z_k_N;
x_k = x_k_N;
v_k = v_k_N;
w_k = w_k_N;


if display_depth ~= 0
  disp('d-step taken ')
end

% PARENTfile: NMS.m