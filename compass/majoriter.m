%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% majoriter: The core of the COMPASS algorithm
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

% Double-loop to solve the MCP:
% 1. minoriter solves linearization (inner loop), yields candidate for new iterate 
% 2. NMS checks acceptance criteria, or 
% 3. finds a new iterate yielding satisfying descent. 
% 4. It is checked if current iterate solves MCP, if not go to 1. (outer loop).  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



take_pg_step = 0;               % flag whether pg step shall be taken next or not
restart = 0;                    % flag whether a restart has to be made  
stop=0;                         % flag for stopping major iteration loop
abort = 0;                      % flag for stopping minor iteration loop  
unsuccessful = 0;               % flag whether pivot procedure was successful or not 

Delta = Delta0;                 % measure of closeness permitting d-steps
d_count = 0;                    % counts d-steps
no_solution_LCP_row = 0;        % counts the unsuccessful pivot procedures in a row.
tgone_count = 0;                % counts how often t>1 in a row in minoriter.
								% prevents majoriter from cycling.  


major=0;                        % counts major iterations   


% solution check:
%%%%%%%%%%%%%%%%%
if Psi_value < target
  stop=1; 
  restart=0;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (stop == 0)
	major=major+1;
  if display_depth ~= 0
	disp(''); disp(''); disp('');  
	major
	disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');  
  end
  if display_depth == 2
	disp('%%%%%%%    START  PIVOT    %%%%%%%%');  
	disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');  
  end

  %%%%%%%%%%%%%%%%%%%%
  minoriter
  % pivot routine solves Approx(x_k)=0 and generates path from x_k to x_k_N (or z_k
  % to z_k_N).  
  %%%%%%%%%%%%%%%%%%%%%
  
  if restart == 1    
	stop = 1;
  else

    %%%%%%%%%%%%%%%%%%%%%%%%%
	NMS
    % accepts z_k_N as next iterate (z_k) or searches for a better one.  
    %%%%%%%%%%%%%%%%%%%%%%%%%%

	Psi_value = Psi(z_k,l,u,f,infty);
	if display_depth ~= 0
	  Psi_value
	end
	f_z_k = f(z_k);   % needed below and in minoriter to determine the linear MCP 

	if d_count == 0   % only if no d-step was taken. 

      % update x, w, and v variables:              	
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  for i = 1:n
		if (z_k(i) <= l(i)) & (f_z_k(i) > 0)
		  x_k(i) = l(i) - f_z_k(i);
		  w_k(i) = f_z_k(i);
		  v_k(i) = 0;
		elseif (z_k(i) >= u(i)) & (f_z_k(i) < 0)
		  x_k(i) = u(i) - f_z_k(i);
		  w_k(i) = 0;
		  v_k(i) = -f_z_k(i);
		else 
		  x_k(i) = z_k(i);
		  w_k(i) = 0;
		  v_k(i) = 0;
		end
	  end

      % update checkpoint bestpoint memory Ref etc.: 
	  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  checkpointz = z_k;
	  checkpointx = x_k; 
	  for i = 1:length(memory)
		if memory(i) == 0
		  memory(i) = infty;
		end
	  end
	  if Psi_value <= min(memory)
		bestpoint = checkpointz;
	  end	  
	  for i = 1:length(memory)
		if memory(i) == infty
		  memory(i) = 0;
		end
	  end
	  [dummy, k] = min(memory);
	  if dummy ~= 0
		[dummy, k] = max(memory);
	  end
	  memory(k) = Psi_value;
	  Ref = max(memory);
	  if display_depth == 2
		memory
	  end

	else    

      % also save checkpoint after d_step:
	  if Psi_value <= max(memory)        
		checkpointz = z_k;
		checkpointx = x_k; 
		mark_check = major;

		for i = 1:length(memory)
		  if memory(i) == 0
			memory(i) = infty;
		  end
		end
		if Psi_value <= min(memory)
		  bestpoint = checkpointz;
		end	  
		for i = 1:length(memory)
		  if memory(i) == infty
			memory(i) = 0;
		  end
		end
		[dummy, k] = min(memory);
		if dummy ~= 0
		  [dummy, k] = max(memory);
		end
		memory(k) = Psi_value;
		Ref = max(memory);
		if display_depth == 2
		  memory
		end
	  end	  
	end  
    

    % save point after checkpoint
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if major == mark_check+1                
	  checkpointzN = z_k;
	  checkpointxN = x_k;
	end

	% abort criteria:
	%%%%%%%%%%%%%%%%%%%
	if major > max_major
	  stop =1;
	  restart = 1;
	  if display_depth ~= 0
		disp('          MAJOR ITERATION LITMIT')
	  end
	end

    % solution check:
	%%%%%%%%%%%%%%%%%
	if Psi_value < target
	  stop=1; 
	  restart=0;
	end


  end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% save statistics for output: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_major(restarts+1)=major;

avg_minor = 0;    
for i = 1:major 
  avg_minor = avg_minor + minor_count(i);
end
if major == 0
  major = 1;
end
save_avg_minor(restarts+1) = avg_minor / major;
save_Psi(restarts+1) = Psi_value;



% stay near last bestpoint in case Psi_value is already considerably small: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Psi_value >= target
if next_x_near == 0
  if (Psi_value < min(max(0.1*Psi_value_0,0.005*(Psi_value_u+Psi_value_l)),0.01*(Psi_value_u+Psi_value_l)))  % is progress already good?
	next_x_near = 1;
	grad_Psi_bestpoint = DPsi(bestpoint,l,u,f,Df,infty)';
	grad_end = bestpoint-grad_Psi_bestpoint; 
    x_next = piB((bestpoint + rand(1,1)*norm(u-l)*0.1* (grad_end-bestpoint)/norm(grad_end-bestpoint)),l,u);
	if display_depth ~= 0
	  disp(' new x0 near last x_k ');
	end
  end
end

if  ((next_x_near == 1) & (restarts >= 3))
  if ((magn(save_Psi(restarts),Psi_value) < 10) & (magn(save_Psi(restarts-1),Psi_value) < 10) & (magn(save_Psi(restarts-2),Psi_value) < 10) )
	x_next = norm(u-l).*rand(1,n); % take a random point as next starting point. 
	if display_depth ~= 0
	  disp('not enough progress:Psi values too close');
	end
  end 
end
end


% PARENTFILE: COMPASSmain.m