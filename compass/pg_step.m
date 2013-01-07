%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pg_step: A projected gradient step may occur in the NMS procedure.   
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

% A pg-step occurs e.g. if the step size was too small in the serach
% during a watchdog step, or if the pivot technique did not yield an
% appropriate point as next iterate.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


grad_Psi_bestpoint = DPsi(bestpoint,l,u,f,Df,infty)';
								% value of gradient of Psi at bestpoint 
s=norm(u-l)/norm(grad_Psi_bestpoint);

Psi_bestpoint=Psi(bestpoint,l,u,f,infty); 
grad_end = bestpoint-grad_Psi_bestpoint; 
								% end point of search direction
search_type_change = 0;         % switch the searchtype
m = 0;
search = 1;


while search == 1  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  alpha_pg = s*(0.4)^m;

  if pg_search_type == 1        % z_alpha can be defined in two ways: 
	z_alpha = (bestpoint + alpha_pg*(piB(grad_end,l,u)-bestpoint));
  elseif pg_search_type == 2  
	z_alpha = piB((bestpoint + alpha_pg*(grad_end-bestpoint)),l,u);
  end

  Psi_z_alpha = Psi(z_alpha,l,u,f,infty);

  if ( grad_Psi_bestpoint*(bestpoint - z_alpha)' < 0 ) 
								% criterion which criterion to apply:	
	    % descent criterion 1
    	if (Psi_z_alpha  >  (Psi_bestpoint - sigma_NMS* grad_Psi_bestpoint*(bestpoint-z_alpha)' ))
		  m = m+1;
		else 
		  search = 0;
		  z_k = z_alpha;
		  if display_depth ~= 0
			disp('% pg-step taken (case 1)')
		  end 
		end
  else
	    % descent criterion 2
    	if (Psi_z_alpha  >  (1-alpha_pg*sigma_NMS)*Psi_bestpoint)
		  m = m+1;
		else 
		  search = 0;
		  z_k = z_alpha;
		  if display_depth ~= 0
			disp('% pg-step taken (case 2)')
		  end 	
		end
  end	      

  if m == max_m_alpha_pg
	if search_type_change == 0  
	  m = 0;
	  if pg_search_type == 1  % change searchtype, try again.    
		pg_search_type = 2;
	  elseif pg_search_type == 2
		pg_search_type = 1;
	  end
	  search_type_change = 1;
	else
	  search = 0;
	  stop = 1;
	  restart=1;
	  if display_depth ~= 0
		disp('alpha_pg  too small in pg_step: Restart. ')
	  end
	  next_x_near =1;
	  x_next = piB((bestpoint + 0.5*rand(1,1)*s*(grad_end-bestpoint)),l,u);
	  if display_depth ~= 0
	    disp('% new x0 near last x_k ');
	  end

	end
  end

end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if search_type_change == 1
  search_type_change = 0;
  if pg_search_type == 1    % change back
	pg_search_type = 2;
  elseif pg_search_type == 2
	pg_search_type = 1;
  end
end


for i = 1:length(memory)
  memory(i) == Psi_z_alpha;
end


%%%%%%%%%%%%%%%%
check_niceness
%%%%%%%%%%%%%%%%


%PARENTFILE: NMS.m
