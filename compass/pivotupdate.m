%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pivotupdate: A new variable and column enters the Basis equation. 
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

% The basis includes t and r at position 'j_old'. 
% The entering variable is called 'Var' of kind 'next', initial value
% 'Var0', covering vector 'Cover' (Var is actually the variable that
% this piece of the path is parametrized by).    

% The value of 'Var' and the position of the new leaving variable 'j'
% are computed. The values of all other variables at the newly obtained
% breakpoint are updated (b_vec, Basis, z_b, w_k, v_k, a_k and t). 
% The variable at position 'j' leaves the basis at its boundary value.  
% t (at its new value) and r are moved from position j_old to j. Var and
% Cover enter the basis at j_old.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


display_min_vec_j =0;
j_old = j;                       % store the current index of t. 
zwvt(j_old) = 0;                 % 0 is the zwvt-value of t (and r)
								 %  r ist hier in der Basis! stelle j_old



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tentative 1st Basis update: Cover enters at j_old: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basis is invertible here!


back_to_start = 1;
while back_to_start == 1

Basis_1 = Basis;

% check if next basis is regular:
d_1 = (Cover'-Basis(:,j_old))'; % difference vector; change in Basis (1st update). 
Basis_1(:,j_old) = Cover;       % Basis_1 includes Cover then. 
delta_1 = (1+H(j_old,:)*d_1');  % if this is not small, the new Basis is regular. 



if abs(delta_1) < 1e-8          % Cannot update H (=inv(B)) properly. 
								% exchange a z-column for an a-column:  

  b_vec_old = b_vec;
  Basis_old = Basis;
  H_old = H;

  %%%%%%%%%%%%
  rankdef
  % returns Basis_new, H_new and b_vec_new; inputs:  _old (respectively)
  % cures rank deficiency. by changing z-columns for a-columns. 
  % updates also other necessary variables like M, M_rem, etc.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Basis = Basis_new;               
  H = H_new;
  b_vec = b_vec_new;

  back_to_start = 1;
else
  back_to_start = 0; 

  % update already now: (as will happen after Cover enters.)
  beta_1 = 1/delta_1;
  H_1 = H-(beta_1.*(H*d_1'))*H(j_old,:);





  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %            take pivot step: 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % compute the coefficient vector, get the value 'Var' of the now
  % entering variable and get the next leaving variable's index 'j':  


  coeff =  (H*Cover')';         % get actual coefficient vector 


  %%%%%%%%%%%%%%%%%%
  check_niceness
  %%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%
  get_min_vec
  % searches for the new value of the entering variable and the index of 
  % the next leaving variable. 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  min_vec_pos = pos(min_vec) + pos(-min_vec);
  [Var,j] = min( min_vec_pos );
  Var = min_vec(j);     


  %%%%%%%%%%%%%%%%%%%%%%%%%
  prevent_move
  % If variables are at bound and would exceed their bounds the procedure is
  % prevented from taking a step. One of those variables leaves instead.
  %%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % check if same z enters and leaves again: 
  % Var, if of kind z, may also hit its own opposite border first, and exit	again
  % just as it enters.      
  if  next == 1
	if ((direction == 1) & (Var > u(j_old)-l(j_old))) 
      Var = u(j_old)-l(j_old);
      j=j_old;    
	end
	if ((direction ==-1) & (Var < l(j_old)-u(j_old))) 
      Var = l(j_old)-u(j_old);
      j=j_old;              
	end          
  end
  %%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ray termination: 
  if Var == ray
	abort = 1; 
    % go a bit on the path (right direction!!), take that point if possible:
	if direction == -1
	  Var = -0.02;       
	elseif direction == 1
	  Var = 0.02;       
	end
	if (b_vec(j_old)-Var*coeff(j_old)) < 0		% new value for t
	  unsuccessful = 1;   	    % no solution to the pivot procedure. 
	  if display_depth == 2
		disp('          ray termination; t value negative along the ray')
	  end
	else 
	  if display_depth == 2
		disp('          ray termination')
	  end
	end						% the updates (rest of file) are made as
							% in any other iteration.    
  end                          
  %%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  % truncate the path if t>1: 
  if (b_vec(j_old)-Var*coeff(j_old)) > 1         % new value for t
	Var = (b_vec(j_old)-1) /coeff(j_old);
	abort =1;
	if display_depth == 2
	  disp('          path truncated at t=1')	
	end
  end 
  %%%%%%%%%%%%%%%



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if abort == 1
	back_to_start = 0;
  else


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%     tentative 2nd Basis update: r enters at j: 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Basis is invertible here!
	% if it is not invertible after this 2nd update I have to redo the
	% first update and go back to the start. 

	% check if next Basis is regular:
	Basis_2 = Basis_1;
	d_2 = (r'-Basis_1(:,j))';     % difference vector; change in Basis (2up). 
	Basis_2(:,j) = r;             % includes Cover (pos j_old) and r (pos j). 
	delta_2 = (1+H_1(j,:)*d_2');  % if not small, Basis is regular after r enters. 

	if abs(delta_2) <1e-8     	% Basis singular after 2nd update. Go back to
								% original Basis and H and try to exchange a
								% z-column for an artificial variable column.  
	  b_vec_old = b_vec;
	  Basis_old = Basis;
	  H_old = H;

	  %%%%%%%%%%%%
	  rankdef
	  % returns Basis_new, H_new and b_vec_new; inputs:  _old (respectively)
      % cures rank deficiency. by changing z-columns for a-columns. 
      % updates also other necessary variables like M, M_rem, etc.
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	  Basis = Basis_new;               
	  H = H_new;
	  b_vec = b_vec_new;

	  back_to_start = 1;
	else 
	  back_to_start = 0;

	  % update already now: (as will happen after r enters.)
	  beta_2 = 1/delta_2;
	  H_2 = H_1-(beta_2.*(H_1*d_2'))*H_1(j,:);
	  
	end
  end 
end 
end 

if display_depth == 2
 	j
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      take the step:
%      update b_vec, Basis, z_b, w_k, v_k, a_k, zwvt and t:  (loop preparation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b_vec = (b_vec-Var*coeff);      % update b_vec. Now: b_vec(j) == 0,l or u (or 1)
Var = Var+Var0; 

if abs((Basis*b_vec')' + Var*Cover  - (-q+r-(M_rem*z_rem')')) > 1e-4 
  abort = 1;
  if display_depth ~= 0
	disp('          ERROR 7') 
  end
end

t = b_vec(j_old);               % update t 
b_vec(j_old) = Var;             % update b_vec (Var replaces t)
zwvt(j_old) = next;             % update zwvt

for i = 1:n
  if zwvt(i) == 1
	z_b(i) = b_vec(i);          % update z_b
  elseif zwvt(i) == 2
	w_k(i) = b_vec(i);          % update w_k
  elseif zwvt(i) == 3
	v_k(i) = b_vec(i);          % update v_k
  elseif zwvt(i) == 4
	a_k(i) = b_vec(i);          % update a_k
  end
end  

% Now t and r are not in the basis. At j_old, the new variable and its
% covering vector entered (Var and Cover). At position j, the variable
% that will now leave the basis, has become u, l, or 0 in b_vec, and in
% the corresponding variable vector (i.e. in z_b, w_k, or v_k).   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     update b_vec, Basis and H:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abort ~= 1 
  % update Basis and H (=inv(Basis)) and b_vec:
  b_vec(j) = t;            
  Basis = Basis_2;          % incl. Cover and r now.  
  H = H_2;
end				   		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The j'th variable just left the basis (but is still left in z_b, w_k
% or v_k).   


%%%%%%%%%%%%%%%%%%
check_niceness
%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t < 0
  if t < -1e-8
	unsuccessful=1;     % if this happens in last pivot iteration:  NO SUCCESS.   
  else
	t = 0;
  	unsuccessful=0;
  end
else
  unsuccessful=0;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

 

% PARENTFILE: pivotloop.m