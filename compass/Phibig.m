%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phibig: Modified Fisher Burmeister function (see Billups, MCP-function)
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

% Finding a zero of this function is equivalent to finding a solution to
% the MCP.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function  Phi  = Phibig(z,l,u,f,infty)

  for i = 1:length(l)
    f_z = f(z);
    if     ((-infty <  l(i)) & (u(i) >= infty))
       Phi(i) = phismall( z(i) - l(i) , f_z(i));

    elseif ((-infty >= l(i)) & (u(i) <  infty))
       Phi(i) = -phismall( u(i) - z(i) , -f_z(i));

    elseif ((-infty <  l(i)) & (u(i) <  infty)) 
       Phi(i) = phismall( z(i) - l(i), phismall( u(i) - z(i), -f_z(i)));

    elseif ((-infty >= l(i)) & (u(i) >= infty))
       Phi(i) = - f_z(i); 

    end
  end
end
