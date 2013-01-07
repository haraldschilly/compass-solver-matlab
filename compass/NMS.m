%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NMS: Non-Monotone Stabilization scheme
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

% descent or distance criteria are used to check whether x_k_N is
% accepted as the new iterate. If it is not, a search is performed in 
% order to determine an acceptable point.  

% possible steps:
% d-steps, m-steps, watchdog steps, projected gradient steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if take_pg_step == 1

  %%%%%%%%%%
  pg_step
  % Returns z_k, new checkpoint
  %%%%%%%%%%%%%%%%%%

  take_pg_step = 0;
  d_count=0;
  Delta = Delta0;
  mark_check = major;

else 

  if (norm(z_k_N - z_k) <= Delta) & (d_count <= max_d_count) % distance criterion

	if (Psi(z_k_N,l,u,f,infty) <= Ref + Ref*0.1) % weak descent criterion

	  %%%%%%%%%
      d_step
      % Returns z_k, x_k, w_k and v_k
	  %%%%%%%%%%%

      d_count = d_count+1;
	  Delta = Delta-deltachange;
	else
 
	  %%%%%%%%%%%%%%%
      watchdog_step
	  % Returns z_k, new checkpoint
	  %%%%%%%%%%%%%%%%%%%

      d_count=0;
	  Delta = Delta0;
	  mark_check = major;
	end

  else

	%%%%%%%%%%%%%
	m_step
	% Returns z_k, new checkpoint
	%%%%%%%%%%%%%%%%%

	d_count=0;
	Delta = Delta0;
	mark_check = major;
  end

end


% PARENTFILE: majoriter.m