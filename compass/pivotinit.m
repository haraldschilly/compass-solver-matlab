%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pivotinit: Initialization of the pivot routine. 
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

% Basis equations are introduced. Rank deficiency is checked. 
% The concept of the pivot procedure is explained. 
% The first step of procedure is taken: t and r enter the basis, at j-th
% position (which is to be determined in get_min_vec). The j-th variable
% leaves the basis at its bound.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


minor=0;          % counts the breakpoints (iterations inside the pivot routine)
				  % first count happens after t enters (pivotinit.m)   
Iv=eye(n);        % These are the matrices used for v and w in the basis
Iw=Iv;				

if abs((M*z_k'-Iw*w_k'+Iv*v_k')' - (-q+r)) > 1e-6
  if display_depth ~= 0
	disp('          ERROR 1') 
  end
  abort = 1;
end							% this must hold at current point


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           construct b_vec and Basis:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct the initial basis matrix and basis vector, called 'Basis' and 'b_vec'
% respectively:  
% As many slack columns as possible are used (unity matrix columns corresponding to w
% or v kind variables with value =0). z kind variables enter the basis only in the
% case when l<z<u holds strictly!  

M_b = zeros(n);                 % Basic part of M (non-basic columns zero)
M_rem = M;                      % remainder of M (basic columns zero)
z_rem = z_k;                    % remainder of z_k (basic variables zero)
zwvt = zeros(1,n);              % vector remembering which index of b_vec is z, w, v,
								% or t.  
for i = 1:n               
  if ((l(i) < z_k(i)) & (z_k(i) < u(i))) 
								% index of kind z: only if   l < z < u      
    zwvt(i) = 1;            

    M_b(:,i) = M(:,i);          % i-th column of M becomes basic
    M_rem(:,i) = zeros(n,1);                
    z_rem(i) = 0;
    Iv(i,i) = 0;                % Make non basic columns of Iv 0 
    Iw(i,i) = 0;                % Make non basic columns of Iw 0 

  elseif z_k(i) == l(i)             % index of kind w also if w is slack (=0)
    zwvt(i) = 2; 
    Iv(i,i) = 0;

  elseif z_k(i) == u(i)             % index of kind v also if v is slack (=0)
	zwvt(i) = 3;
    Iw(i,i) = 0;     
  elseif ((z_k(i) < l(i)) | (u(i) < z_k(i))) 
	if display_depth ~= 0
	  disp('followup ERROR pivotinit:    z_k(i)  not in box: ');
	end
  end
end

z_b = z_k - z_rem;              % basic part of z_k (non-basic variables zero) 
a_k = zeros(1,n);               % artificial variable vector:
Ia = zeros(n);                  % covering vectors for artivficial variables. 
b_vec = z_b + w_k + v_k + a_k;  % vector of all basic variables (ordered!)
Basis = M_b - Iw + Iv + Ia;     % Matrix of all basic columns (ordered!)
 
if abs((Basis*b_vec')' - (-q+r-(M_rem*z_rem')')) > 1e-6
  if display_depth ~= 0
	disp('          ERROR 3') 
  end	
end  






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Incase Basis is rank deficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if abs(rcond(Basis)) < 1e-8
  if display_depth ~= 0
	disp(' Initial Basis rankdef.');
  end

  % Beginning from the unity matrix, it is checked, which columns of
  % Basis may enter so that the resulting matrix (then taken to be the
  % basis) is invertible. All columns that can not enter this new basis
  % will be replaced by artificial unity matrix columns and
  % corresponding variables will be replaced by artificial variables of
  % kind a:

  Basis_1 = eye(n);
  H_1 = Basis_1;  
  is_z = [];
  d_ = zeros(n);
  for i = 1:n
	if zwvt(i)==1
	  is_z = [is_z, i];         % vector saving the z-indices. 
	elseif zwvt(i)==2
	  Basis_1(i,i)=-1; % w-indices are negative in Basis! 
      H_1(i,i)=  -1;
	end 
  end
  II=eye(n);
  d_(:,is_z) = Basis(:,is_z)-II(:,is_z); 
								% difference vectors (replace as many as possible!)  

  complete = 0;
  while complete == 0

	delta = zeros(1,n);
	for i = 1:length(is_z)
		delta(is_z(i)) = (1+H_1(is_z(i),:)*d_(:,is_z(i))); 
								% denominator of each possible z-variable
	end

	[delta_enter, j_enter] = max(abs(delta));
	if delta_enter < 1e-8
      complete = 1;
	else
	  delta_enter=delta(j_enter); %save with right sign!
	  Basis_1(:,j_enter)=Basis(:,j_enter); % update Basis
	  beta = 1/delta_enter;
	  H_1 = H_1-beta.* H_1*d_(:,j_enter)*H_1(j_enter,:); % update inverse basis.
	  if length(is_z) == 1
		complete = 1;
		if display_depth == 2
		  disp('no a inserted in cure');
		end
	  end

	  if length(is_z) > 0
		kill_j_enter=1;
		k = 1;
		while kill_j_enter == 1
		  if is_z(k) == j_enter
			is_z(k)=[];		  
			kill_j_enter=0;
		  else
			k = k+1;
		  end
		end
	  end
	end
  end

  zwvt(is_z) = 4;                   % remaining z-indices become a's:
  Ia(is_z,is_z)=1; 
  M_rem(:,is_z) = M(:,is_z);
  z_rem(is_z) = z_b(is_z);          % now in M_rem but not at bound! 
  M_b(:,is_z) = zeros(n,is_z);                  
  z_b(is_z) = 0;

  if display_depth == 2
	disp('artificials inserted at position:');
	art_pos = is_z
  end

  b_vec = z_b + w_k + v_k + a_k;
  Basis = M_b - Iw + Iv + Ia;   
  H = H_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



else 
  H = inv(Basis);
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Start of the pivot procedure:    t*r enters the basis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


back_to_start = 1;
while back_to_start == 1

  coeff = (H*r')';              % coefficients determine how the variables (b_vec(i)) 
					 			% change as t increases.          


  % abs((Basis*b_vec')' = (-q+r-(M_rem*z_rem')'))  ... must hold for
  % small t. Now t will increase until one index of 

  % b_vec - t*coeff
  
  % becomes zero. This index will be called 'j' in the following. The
  % '=' sign holds during the whole procedure.   


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % (   Short explanation of the further procedure:

  % Later b_vec will be updated as    
  %     b_vec = b_vec - t*coeff, 
  % r will enter the 'Basis' matrix instead of the j-th column, and t
  % will enter the 'b_vec' vector instead of the j-th variable. Incase
  % the j-th variable was of kind 'z', it will be moved to the
  % 'M_rem*z_rem' part of the equation (at its new value, i.e. its
  % boundary), if it was of kind 'a', 'w' or 'v' it will simply
  % disappear (boundary values of 'w' or 'v' kind variables are zero).   

  % According to the 'pivot rules', another kind of variable will enter
  % the basis in the next step (pivotloop.m), just as t*r enters here
  % (this file), and cause another variable (most probably at a
  % different index position) to leave.    

  % As a convention, variables will never change their index positions.
  % Instead, t and r will jump to the position that just became empty in
  % each iteration step.   

  % The algorithm terminates successfully if t leaves the basis at its boundary t=1.  
  % ) 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  direction = 1;

  %%%%%%%%%%%%%%%%%
  check_niceness
  %%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%
  get_min_vec
  % searches for the new value of the entering variable and the index of 
  % the next leaving variable. 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  [t,j] = min(pos(min_vec)+pos(-min_vec));
  t = min_vec(j);     
								% value of t and the index j, in which the next Basis
								% and variable update will occur.   
								% Choose t over the positive values of the min_vec
								% components, but save t with the right sign!

  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  prevent_move
  % If variables are at bound and would exceed their bounds the procedure is
  % prevented from taking a step. One of those variables leaves instead.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
   

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if t >= 1
	abort = 1; 
	t=1; 
  	tgone_count = tgone_count+1; % prevents the outer algorithm from cycling. 
  else
  	tgone_count = 0;
  end
  if tgone_count >= 3;
  	tgone_count = 0;
  	unsuccessful = 1;
  	take_pg_step = 1;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%
  if abs((Basis*(b_vec-t*coeff)')'+t*r - (-q+r-(M_rem*z_rem')')) > 1e-4
	abort = 1;
	if display_depth ~= 0
	  disp('          ERROR 6') 
	end
  end
  %%%%%%%%%%%%%%%%



  if abort == 1
	back_to_start = 0;
  else

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tentative update: r enters at pos.j
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Basis is invertible here!

    % check if next Basis is regular:
	Basis_2 = Basis;            % tentative basis after the update  
	d_2 = (r'-Basis_2(:,j))';   % difference vector; change in Basis. 
	Basis_2(:,j) = r';
	delta_2 = (1+H(j,:)*d_2');    

	if abs(delta_2) <1e-8
      % Basis would be singular after this update!

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
    
      % update H_2 already now: (as will happen after r enters)
	  beta_2 = 1/delta_2;
	  H_2 = H-beta_2.*H*d_2'*H(j,:);
	end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           take the step:
%           update b_vec, Basis, z_b, w_k, v_k:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
							
Var = t;                        % used in the loop like this 						
Cover = r;                      % used in the loop like this
b_vec = (b_vec-Var*coeff);      % update b_vec

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

% The corresponding point (z_b, w_k, v_k, a_k) is the first breakpoint
% of the path. The variable that will now leave the basis, has become u,
% l, or 0 in b_vec, and in its corresponding variable vector (i.e. in
% z_b, w_k, or v_k). 

if abs((Basis*(b_vec)')'+t*r  - (-q+r-(M_rem*z_rem')')) > 1e-4
  if display_depth ~= 0
	disp('          ERROR 6a') 
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     update b_vec, Basis and H:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abort ~= 1 

  % update Basis and H (=inv(Basis)) and b_vec:
  b_vec(j) = t;
  Basis = Basis_2;
  H = H_2; 						% The j'th variable just left the basis (but is still
								% left in z_b, w_k or v_k).  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if display_depth == 2
 	j
end


%%%%%%%%%%%%%%%%%%
check_niceness
%%%%%%%%%%%%%%%%



% PARENTFILE: minoriter