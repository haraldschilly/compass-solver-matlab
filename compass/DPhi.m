%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DPhi: Jacobian of Phibig
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

function  DPhi = DPhi(z,l,u,f,Df,infty)
  f_z = f(z); 
  Df_z=Df(z);
  for i = 1:length(l)
	for k = 1:length(l)

        
      if     ((-infty <  l(i)) & (u(i) >= infty))
		if i==k
		  DPhi(i,k) = (d_1_phi(z(i)-l(i),f_z(i)) +  d_2_phi(z(i)-l(i),f_z(i))*Df_z(i,k)); 
		else
		  DPhi(i,k) = d_2_phi(z(i)-l(i),f_z(i))*Df_z(i,k);
		end


      elseif ((-infty >= l(i)) & (u(i) <  infty))
		if i==k
		  DPhi(i,k) = (d_1_phi( u(i) - z(i) , -f_z(i))  + d_2_phi( u(i) - z(i) ,  -f_z(i))*Df_z(i,k));			  
		else
		  DPhi(i,k) = d_2_phi( u(i) - z(i) , -f_z(i))*Df_z(i,k);
		end


      elseif ((-infty <  l(i)) & (u(i) <  infty)) 
		if i==k
		  DPhi(i,k) = (d_1_phi( z(i) - l(i), phismall( u(i) - z(i), -f_z(i))) -  (d_2_phi( z(i) - l(i), phismall( u(i) - z(i), -f_z(i))) *d_1_phi( z(i) - l(i), phismall( u(i) - z(i), -f_z(i)))) - (d_2_phi( z(i) - l(i), phismall( u(i) - z(i), -f_z(i)))*Df_z(i,k))^2);
		else
		  DPhi(i,k) = (-(d_2_phi( z(i) - l(i), phismall( u(i) - z(i), -f_z(i)))*Df_z(i,k))^2);
		end


      elseif ((-infty >= l(i)) & (u(i) >= infty))
		-Df_z(i,k);
      end


	end
  end

end