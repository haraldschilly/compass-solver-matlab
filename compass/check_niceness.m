%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check_niceness: routine controlling if variables are in bounds. 
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

% This routine checks if all variables are nice (i.e. inside their
% bounds); according to the theory all indices should always be inside
% their bounds. Rounding errors could lead to a numerical violation of
% the bounds. If any variables' values are out of their bounds they are
% projected on their bounds in this routine.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n
  if in_pivot == 1      % file is called in the pivot routine

    %check only z_b, z_rem is ok! 
	if zwvt(i) == 1  
	  if z_b(i) < l(i)          % Ignore nonbasic z's 
		z_b(i) = l(i);
	  elseif u(i) <  z_b(i)
		z_b(i) = u(i);
	  end
	end 
 
  else                  % file is called outside the pivot routine

	%check z_k
	if zwvt(i) == 1  
	  if z_k(i) < l(i)          % Ignore nonbasic z's 
		z_k(i) = l(i);
	  elseif u(i) <  z_k(i)
		z_k(i) = u(i);
	  end
	end  
  end

  if w_k(i) < 0
	w_k(i) = 0;
  end
  if v_k(i) < 0
	v_k(i) = 0;
  end  
  if a_k(i) < 0
	a_k(i) = 0;
  end

end


% Parentfiles: pivotupdate.m, pivotinit.m, minoriter.m, NMS.m
