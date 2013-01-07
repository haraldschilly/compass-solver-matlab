%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crash: The crash technique used in COMPASS
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

% crash determines a good approximation to the the active set at the solution. 
% This routine starts at x0, takes a few projected Newton steps, returns z_k.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



no_change_row = 0;
crash_count = 0;                  % counts the iterations in crash procedure  
if n <= min_dim_crash             % minimum problem dimension to perform crash
  do_crash = 0;
end

if do_crash == 1
  if display_depth ~= 0
	disp('enter crash phase');
	disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %             The Reduced system: 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

								% The index vector Red (below) contains all inactive
								% indices, and those active 
								% indices with wrong value of f.
								% the search (later) is restricted to these indices,
								% i.e. to the face of B defined by 
								% Red.  

								% starting from z_k, the aim is to find d_k (Newton
								% direction of reduced system)  
								% and search for alpha so that
								%        z_k+1 := piB(z_k - alpha*d_k)    
								%
								% yields descent in the NORMALMAP. z_k is returned
								% after at most 'max_crash_count' steps, if the  
								% stepsize in the search is too small, or if not
								% enough indices in the active set change.  

  Ind = (1:n);                  % define index vector Ind containing all indices 1 to
								% n.  

  fixed_ind = zeros(n,1);       % index vector with 1 entries at index positions
								% corresponding to fixed indices (active with right f
								% value), and 0  entries at index position
								% corresponding to indices in the face defined by
								% Red.   
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           The Crash loop: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


while (do_crash == 1)
  crash_count = crash_count+1;

  Red = Ind;                    % define Red 
  fixed_ind_old = fixed_ind;
  
  f_z_k = f(z_k);            	% f at z_k     
  z_k_red = z_k;             	% define reduced z_k
  f_red_z_k = f_z_k;         	% define reduced f at z_k
  Jac_f_red_z_k = Df(z_k);		% define reduced Jacobian matrix of f at x_k 

  % reduce:
  for i = 1:n                            
	if (((l(n+1-i)==z_k(n+1-i)) & (f_z_k(n+1-i)>=0)) ||	((u(n+1-i)==z_k(n+1-i)) & (f_z_k(n+1-i)<=0)))
	  Red(n+1-i)=[]; 
	  fixed_ind(n+1-i) = 1;
	  z_k_red(n+1-i) = [];           
	  f_red_z_k(n+1-i) = [];
	  Jac_f_red_z_k(:,n+1-i) = [];
	  Jac_f_red_z_k(n+1-i,:) = [];
	end
  end


  if display_depth ~= 0
	active_change = (norm(fixed_ind_old - fixed_ind))^2
  else
	active_change = (norm(fixed_ind_old - fixed_ind))^2;
  end							% number of index changes in the active set (active
								% with right f value) 


  % In case of rank deficiency: 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ((rcond(Jac_f_red_z_k) < 1e-16)  )
	epsilon = 1e-16;            % a perturbation parameter is used to make the matrix
								% numerically nonsingular (size of Epsilon?):
	while ((rcond(Jac_f_red_z_k) < 1e-16)  )
	  epsilon = epsilon*10;
	  Jac_f_red_z_k = Jac_f_red_z_k + epsilon*eye(length(Red));
	end
	if display_depth ~= 0
	  disp('reduced grad_f(z) rank deficient in crash. Perturbation parameter:');
	  epsilon
	end
  end



  % derive search direction d_k for the search: 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  d_k_red = Jac_f_red_z_k \ f_red_z_k'; % Newton direction of reduced system
  d_k = zeros(1,n);             
  for i = 1:length(Red)
	d_k(Red(i)) = d_k_red(i);	% search direction d_k of non reduced system. 
  end





  % search for alpha:
  %%%%%%%%%%%%%%%%%%%

  searchalpha = 1;              % flag for doing the search for alpha. 
  m = 0;                                             

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while (searchalpha == 1)
  								% z_k is given, but the search for alpha will be done
								% in terms of x_k.  
								% Values for w and v are obtained by projecting -f(z)
								% onto the normal cone of B at z.  
								% Then x_k:=z_k-w+v, and z_k = piB(x_k): 

	alpha = (0.5)^m;
	z_alpha = piB(z_k - alpha*d_k, l, u);      						
	f_z_alpha = f(z_alpha);
	w_alpha = zeros(1,n);         
    v_alpha = w_alpha; 

	% assign v and w the projection of f(z) on the respective index: 
	for i = 1:n
	  if ((z_alpha(i) == u(i)) & (f_z_alpha(i) <= 0))
		v_alpha(i) = -f_z_alpha(i);
	  end	
	  if ((z_alpha(i) == l(i)) & (f_z_alpha(i) >= 0))
		w_alpha(i) = f_z_alpha(i);
	  end
	end 

	x_alpha = z_alpha - w_alpha + v_alpha;
	% x_alpha is needed for the search: 

	if (norm(normalmap(x_alpha,l,u,f)) < (1-sigma_crash*alpha)* norm(normalmap(x_k,l,u,f))) 
	  searchalpha = 0;
	else
	  m = m+1;
	  if m == max_m_alpha_crash      % check if stepsize becomes too small 
		searchalpha = 0;
		do_crash = 0;
		if display_depth ~= 0
		  disp('          Step size too small in alpha-search')
		end
	  end
	end

  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % new iterate found:
  %%%%%%%%%%%%%%%%%%%%

  z_k = z_alpha;
  x_k = x_alpha;


  % abort criteria: 
  %%%%%%%%%%%%%%%%%

  if (active_change < min_active_change)
	no_change_row = no_change_row+1;
	if no_change_row > max_no_change_row
	  do_crash = 0;
	  if display_depth ~= 0
		disp('          too many steps with insufficient change in active set.')
	  end
	end
  else
	no_change_row=0;

  end

  if crash_count == max_crash_count
	do_crash = 0;
	if display_depth ~= 0
	  disp('          Maximal possible steps.')
	end
  end

  if do_crash == 0
	if display_depth ~= 0
	  disp('          Crash procedure terminated. Number of steps taken:')
	  crash_count
	end
  end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% end crash loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Psi_value = Psi(z_k,l,u,f,infty);
if display_depth ~= 0
  Psi_value
end


% accept crash result if it is not too much higher than Psi_value at x0: 
if Psi_value > Psi_value_0 * 1.2
  x_k = x0;
  z_k = piB(x_k,l,u);
  Psi_value = Psi_value_0;
  if display_depth ~= 0
	disp('x0 used');
  end
else
  if display_depth ~= 0
	disp('crash accepted');
  end
end


% update x, w, and v variables for use in majoriter:
f_z_k = f(z_k);
for i = 1:n
  if (z_k(i) <= l(i)) & (f_z_k(i) > 0)  % '=' never occurs due to 'piB'.
	x_k(i) = l(i) - f_z_k(i);
	w_k(i) = f_z_k(i);
	v_k(i) = 0;
  elseif (z_k(i) >= u(i)) & (f_z_k(i) < 0) % neither here. always at z_k. 
	x_k(i) = u(i) - f_z_k(i);
	w_k(i) = 0;
	v_k(i) = -f_z_k(i);
  else 
	x_k(i) = z_k(i);
	w_k(i) = 0;
	v_k(i) = 0;
  end
end
a_k = zeros(1,n);



% save checkpoints: 
checkpointx = x_k; 
checkpointz = z_k;
bestpoint = checkpointz;
mark_check = 0;                 % remembers iteration number of last checkpoint

Ref = Psi_value;                % Reference value
memory = [Ref 0 0 ];            % remembers best last three Psi values. 




%PARENTfile: COMPASSmain.m