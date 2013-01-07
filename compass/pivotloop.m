%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pivotloop: The loop (inner loop) and the pivot rules.  
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

% The kind of the j-th variable is determined (just left the basis).  
% According to specific pivot rules the entering variable and its
% covering vector are determined.  

% M_rem z_rem, z_b are updated completely. 
% The pivot routine is called for the next step (pivotupdate).  
% A solution to the linear MCP is returned. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cycling = 0;                    % flag for cycling
cycle_det = zeros(n,1);         % cycling detection: counts how many times (value)
								% a variable change happens at the same position 
remember_t =[];                 % stack remembering t values in case of cycling.  


while abort == 0                % loop starts   
  minor = minor+1;

  switch zwvt(j)                % the j-th variable just left the basis. 
	case 1                      % the j-th variable was of kind 'z'

	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  if abs(z_b(j) - u(j)) < 1e-5
 								% z_j left at u_j 

		next = 3;               % v_j enters next at value 0 
		cycle_det(j) = cycle_det(j)+1;           
		Var0 = 0;
        II=eye(n);    
		Cover = (II(j,:));      % its covering vector 
		direction = 1;          % the variable must rise  
		z_b(j) = 0;             % z_b is updated now completely.
		M_rem(:,j) = (M(:,j));  % update M_rem (needed for ERROR check)
		z_rem(j) = u(j);        % update z_rem 

		pivotupdate

	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  elseif abs(z_b(j) - l(j)) < 1e-5
								% z_j left at l_j 

		next = 2;               % w_j enters next at 0
		cycle_det(j) = cycle_det(j)+1;           
		Var0 = 0;
        II=eye(n);
		Cover = (II(j,:));      % its covering vector 
		direction = 1;          % the variable must rise  
		z_b(j) = 0;             % z_b is updated now completely.
		M_rem(:,j) = (M(:,j));  % update M_rem (needed for ERROR check)
		z_rem(j) = l(j);        % update z_rem 

		pivotupdate

	  else
		if display_depth ~= 0
		  disp('         ERROR.....: z left but not at bound ')
		end
		abort=1;
	  end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	case 2                      % the j-th variable was of kind 'w'

	  next = 1;                 % z_j enters next at lower bound
	  cycle_det(j) = cycle_det(j)+1;           
	  Var0 = l(j);
	  Cover = (M(:,j))';        % its covering vector 
	  direction = 1;            % the variable must rise  
	  M_rem(:,j) = zeros(n,1);  % update M_rem (needed for ERROR check)
	  z_rem(j) = 0;             % update z_rem 
	  
	  pivotupdate

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	case 3                      % the j-th variable was of kind 'v'

	  next = 1;                 % z_j enters next at upper bound
	  cycle_det(j) = cycle_det(j)+1;           
	  Var0 = u(j);
	  Cover = (M(:,j))';        % its covering vector 
	  direction = -1;           % the variable must decrease
	  M_rem(:,j) = zeros(n,1);  % update M_rem (needed for ERROR check)
	  z_rem(j) = 0;             % update z_rem 

	  pivotupdate

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	case 4                      % the j-th variable was of kind 'a'

	  next = 1;                 % z_j enters next at its current value
	  cycle_det(j) = cycle_det(j)+1;           
	  Var0 = z_rem(j);                
	  Cover = (M_rem(:,j))';    % its covering vector 
	  direction = 0;            % the direction of the entering variable can be
								% chosen. 
	  M_rem(:,j) = zeros(n,1);  % update M_rem (needed for ERROR check)
	  z_rem(j) = 0;             % update z_rem 

	  pivotupdate

  end

  % If an 'a' kind variable leaves, it disappears too, and the
  % respective z kind enters; in this case, since this z is not at a
  % bound, the direction can be chosen to increase t.   


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %             abort conditions: 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if abs(t-1) < 1e-8
	abort = 1;                   
  end

  if minor > max_minor
	abort = 1;
	if display_depth == 2
	  disp('          MINOR ITERATION LIMIT')
	end
	unsuccessful = 1;
  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % determine if the procedure cycles:

  if ((minor > maybe_cycle) & (cycling == 0) & (abort == 0))
	if max(cycle_det) > definitive_cycling           
      % cycle_det keeps number of how many times each index enters the basis  

	  cycling = 1;
	  m_cycle = 0;              % counts the next 5 iterations before the maximum t
								% value is chosen.  
	end
  end

  % react to cycling: 
  
  % The point along the cycle with the largest t value
  % is chosen to be the next iterate. 5 more iterations
  % are taken, of which the maximum t value is chosen.
  % (It is assumed that not more than 5 variables are
  % in a cycle.) Minor iterations go on until t reaches
  % that maximum t value again. Then it stops. If all t
  % values are negative, no usable solution can be
  % found. -> pg_step or restart 

  if cycling == 1 

	if m_cycle < cyclesize
	  m_cycle = m_cycle+1;
	  remember_t(m_cycle) = t; 

	elseif m_cycle == cyclesize
	  m_cycle = m_cycle+1;
	  remember_t(m_cycle) = t;
	  t_max = max(remember_t); 
	  
	  if t_max < 0
		abort = 1; 
		if display_depth == 2
		  disp('          cycling; maximum t value negative');
		end
		unsuccessful = 1;		% no solution to the pivot procedure. 
	  elseif  abs(t_max) < 1e-5
		abort = 1; 
		if display_depth == 2
		  disp('          cycling; maximum t value zero');
		end
		unsuccessful = 1;		% no solution to the pivot procedure. 
	  end


	elseif m_cycle > cyclesize
	  if abs(t-t_max) < 1e-4
		if display_depth == 2
		  disp('          cycling; path truncated at maximum t value');
		end
		abort = 1; 				% The current   z_b + z_rem   is returned to NMS. 
		if t <= 1e-5
		  if display_depth == 2
			disp('t = 0');
		  end
		  unsuccessful = 1;
		end
	  end
	end

  end
end


% PARENTFILE: minoriter.m