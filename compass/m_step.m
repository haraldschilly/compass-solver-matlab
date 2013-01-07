%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m_step: A step kind in the NMS procedure in COMPASS.  
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

% Is taken if the distance criterion for a d-step is not satisfied 
% or if d_count is above its limit. The descent criterion is checked. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


grad_Psi_z_k = DPsi(z_k,l,u,f,Df,infty)';
								% value of gradient of Psi at z_k

if ((grad_Psi_z_k * (z_k' - z_k_N')) < 0) 
								% the case where we walk in direction of the gradient

    if Psi(z_k_N,l,u,f,infty)  <=  	Ref - sigma_NMS* grad_Psi_z_k*(z_k - z_k_N)'
	  z_k = z_k_N;
	  if display_depth ~= 0
		disp('m-step taken (case_1)')
	  end
	else 

	  %%%%%%%%%%%%%%%
	  watchdog_step
      % Returns z_k, new checkpoint
	  %%%%%%%%%%%%%%%%%%%

	end

else

    if Psi(z_k_N,l,u,f,infty)  <= 	(1-sigma_NMS)* Ref	
	  z_k = z_k_N;
	  if display_depth ~= 0
		disp('m-step taken (case_2)')
	  end
	else 

	  %%%%%%%%%%%%%%%
	  watchdog_step
      % Returns z_k, new checkpoint
	  %%%%%%%%%%%%%%%%%%%

	end

end 

% PARENTFILE: NMS.m
								









