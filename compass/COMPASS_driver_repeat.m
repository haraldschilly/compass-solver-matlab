%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPASS_driver_repeat: calls COMPASS multiple times on random problems.   
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

% A problem class can be specified in COMPASS_driver_no_clear.m;
% problems of that kind should be randomly generated in that file. 

% 100 problems will be generated and solved. The success rate, as well
% as the average time per successful and unsuccessful problem are
% displayed as results. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


clear;


% how many problems should be solved?
nn = 100;                    

seconds = zeros(1,nn);
av_time_good = 0;
av_time_bad = 0;

ok = zeros(1,nn);
COMPASSdrivecount = 0;


while COMPASSdrivecount < nn
  COMPASSdrivecount = COMPASSdrivecount+1


  %%%%%%%%%%%
  COMPASS_driver_no_clear
  %%%%%%%%%%%%


  seconds(COMPASSdrivecount) =  duration_in_seconds;
  if Psi_value < target
	ok(COMPASSdrivecount) = 1;
	av_time_good = av_time_good + duration_in_seconds;
  else
	av_time_bad = av_time_bad + duration_in_seconds;
  end

end

good= sum(ok);
success_rate = good/nn
av_time_good = av_time_good/good 
av_time_bad = av_time_bad/(nn-good) 






