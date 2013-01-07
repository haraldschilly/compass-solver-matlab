%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prevent_move: Prevents pivot steps in which variables would violate bounds 
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

% If some variables are at a bound and would exceed their bounds, the
% procedure does not move, Var stays == 0. One of these variables (at
% bound that would move beyond it) is chosen to leave instead (index j):
% if neither t nor an a is possible to leave, any variable may leave.  

% variables to leave in this case are chosen in a perpetuating manner in
% order to avoid cycling. Leaving variables when no move is made are
% remembered, new leaving variables are chosen each time if possible.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if move == 0
  Var = 0;
  j = xx_alternative(1);

  t_possible = 0;               % flag whether t may leave or not. 
  a_possible = 0;               % flag whether an a kind may leave or not. 

  for kk = 1:length(xx_alternative)
	kkk = length(xx_alternative)+1-kk;   
								%run through for loop in reverse order.  
	if zwvt(xx_alternative(kkk)) == 0
	  j = xx_alternative(kkk);    
	  t_possible = 1;
	end
	if ((t_possible == 0) & (zwvt(xx_alternative(kkk)) == 4))
	  a_possible = 1;
	  j = xx_alternative(kkk);    % take the smallest a-kind j possible. 
	end
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ((t_possible == 0) & (a_possible == 0))
	no_move_count = no_move_count+1; 
  	if no_move_count > 1
	  k=1;
	  j_already_used = 1;       % flag if j was already used. 

	  while j_already_used == 1
		j_already_used = 0;
		for kk = 1:length(remember_j)         
		  if j == remember_j(kk)% stack, remembers already used j's. 
			j_already_used = 1;
		  end
		end
		if ((j_already_used == 1) & (k < length(j_alternative)))
		  k=k+1;
		  j = j_alternative(k); % take the smallest possible new j to leave. If they
								% were all already used, take the first and, clear
								% the memory.
		elseif ((j == j_alternative(length(j_alternative))) & (j_already_used == 1))
		  remember_j = [];
		  j = j_alternative(1);
		else 
		  j_already_used = 0;
		end
	  end 

	end
	remember_j = [remember_j j];
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ((t_possible == 0) & (a_possible == 1))
	no_move_a_count = no_move_a_count+1; 
	if no_move_a_count > 1
	  k=1;
	  a_already_used = 1;       % flag if this a was already used. 

	  while a_already_used == 1
		a_already_used = 0;
		for kk = 1:length(remember_a)         
		  if j == remember_a(kk)% stack, remembers already used a's. 
			a_already_used = 1;
		  end
		end
		if ((a_already_used == 1) & (k < length(a_alternative)))
		  k=k+1;
		  j = a_alternative(k); % take the smallest possible new a to leave. If they
								% were all already used, take the first and, clear
								% the memory.
		elseif ((j == a_alternative(length(a_alternative))) & (a_already_used == 1))
		  j = a_alternative(1);
		  remember_a = [j];
		else 
		  a_already_used = 0;
		end
	  end 

	end 
	remember_a = [remember_a j];
  end 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if minor == 0
	t = Var;
  end
end

%PARENTFILE: pivotinit.m and pivotupdate.m