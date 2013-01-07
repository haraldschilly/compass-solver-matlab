%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minoriter: The pivot routine of COMPASS
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

% from the current point z_k, a path is constructed with the pivot
% technique and the Newton point z_k_N is obtained as its end point. 
 
% This holds at the current point x_k:   
%   l  =<  z_k  =<  u
%   0  =<  w_k 
%   0  =<  v_k 
%   x_k  ==  z_k - w_k + v_k   (comprised)     
%   v_k*w_k'  ==  0     (complementary)   

%  the path p(t) from  x_k  to  x_k_N  must satisfy: 
%      x_k(t)  ==  p(t)  ==  z_k(t) - w_k(t) + v_k(t)     
							
% The linear MCP is defined by M and q:  find z,w,v, such that
%   M*z - I*w + I*v  ==  -q 
% where at the current point 
%   M*z_k - I*w_k + I*v_k  ==  -q + r    holds. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



no_move_count = 0;              % counts how many times a move is prohibited.  
remember_j = [];                % remembers leaving variables when no move is made. 
no_move_a_count = 0;            % as above with an 'a' leaving.  
remember_a = [];                % remembers leaving 'a's when no move is made. 
up_down = 1;                    % direction switch whenever it is not determined. 

in_minoriter = 1;
in_pivot = 1;					% flag: one is in the pivot procedure files 
			 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    derive the linear MCP:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%
% compute Matrix M:
M = Df(z_k);

% In case of rank deficiency: 
if ((rcond(M) < 1e-16)   )
  epsilon = 1e-16; 				% perturbation parameter as in crash.m				
  while ((rcond(M) < 1e-16)   )
	epsilon = epsilon*10;
	M = M + epsilon*eye(n);
  end
  if display_depth == 2
	disp(' ');
	disp('LCP rank deficient. Perturbation parameter: ');
	epsilon
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% vector q and residual r:
q = f_z_k - (M*z_k')';          
r = normalmap(x_k,l,u,f);   


if Approx(x_k,M,q,l,u) ~=  normalmap(x_k,l,u,f)	% must be equal at current point
  if display_depth ~= 0
	disp('          ERROR  0') 
  end
end





%%%%%%%%%%%%%%%%%
pivotinit
% Approximation is computed, linear MCP is stated, first step in the pivot routine: 
% t and r enter the basis.  
%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%
pivotloop
% the pivot routine, according to the pivot rules. 
%%%%%%%%%%%%%%%%


abort=0;
in_pivot = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solutions of the pivot routine (Newton point): 
z_k_N = z_b + z_rem;                   
w_k_N = w_k;
v_k_N = v_k;
x_k_N = z_k_N - w_k_N + v_k_N;



%%%%%%%%%%%%%%%
check_niceness
%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save point after checkpoint:
								% (necessary if a search must be made directly after
								% a checkpoint was saved)   
if major == mark_check+1
  checkpointzN = z_k_N;
  checkpointxN = x_k_N;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display statements / pivot analysis: 

minor_count(major) = minor;     % saves number of minor iterations for each major
								% iteration  
if unsuccessful == 0
  no_solution_LCP_row = 0;
  if display_depth == 2
	t
	disp('% pivot SUCCESSFUL ')
  end

elseif unsuccessful == 1
  no_solution_LCP_row = no_solution_LCP_row +1;
  if display_depth == 2
	disp('% NO SUCCESS in pivot ')
  end                 		% go back to last checkpoint and take pg step; 
    							% restart if this happens more than
								% max_no_solution_LCP_row times in a row:   
  if no_solution_LCP_row > max_no_solution_LCP_row 
	restart = 1; 
	stop = 1;
	if display_depth == 2
	  disp('no solution to LCP more than 5 times in a row: Restart.')
	end
  else
	take_pg_step = 1;
  end 
  unsuccessful = 0;
end
if display_depth == 2
  disp(' ');
end
in_minoriter = 0;



% PARENTFILE: majoriter.m
