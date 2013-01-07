%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_min_vec: routine computes new entering variables index and value.
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

% search for the new value of the entering variable and the index of the
% new leaving variable. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if direction == 0               
  if  coeff(j_old) < 0          % a left, entering z is not at bound: decide whether
								% z shall rise or fall, so that t rises.  
	direction = 1;              % the update later looks like this:  
								% b_vec = (b_vec-Var*coeff) 
  elseif coeff(j_old) > 0
	direction = -1;
  else 
	if up_down == 1
	  direction =1;
	  up_down =0;
	else
	  direction =-1;
	  up_down = 1;
	end
  end
end

% the role of the variable 'direction' in this for loop: 
% Var starts at zero and has to go into a given direction (rise or decrease:
% direction = 1  or  -1 ). If the value of min_vec(i) has the wrong sign (depending
% on e.g. coeff(i) ), it shall not be chosen to be the next index for a leaving
% variable. Hence the 'ray' value assignment.  


min_vec = zeros(1,n);
ray = max(max(abs(u)), max(abs(l)))*1000;       
j_alternative = [];
a_alternative = [];
xx_alternative = [];
move = 1;

% min_vec is the vector that contains for each index 'i' the weighted (coeff(i))
% distance  of the i-th basic valriable's value and its boundary, i.e. the value that
% Var must attain so that the i-th variable hits its bound.
 
% The minimum element of min_vec is hence the new value of 'Var'. 

% min_vec is constructed in the following 'for' statement. If a variable, say at
% position 'k', will increase (or decrease) infinitely as Var increases, min_vec(k)
% will be set to a very large value (=ray).   


% display purpose:
a_zero = [];
ds_v=0;
V_ZERO_and_FALLS = [];
ds_w=0;
W_ZERO_and_FALLS = [];
ds_a=0;
A_ZERO_and_FALLS = [];
ds_zu=0;
Z_U_and_RISES = [];
ds_zl=0;
Z_L_and_FALLS = [];




% calculate the (weighted) 'distance' to bound of...  




for i = 1:n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (zwvt(i) == 1)             % ...each z variables in the basis
    if coeff(i) ~= 0

%%%%%%%%%%%%%%%%    lower bound    %%%%%%%%%%%%%%%%%%%%%%%%        

	  if (((direction == 1) & (0 < coeff(i))) | ((direction == -1) & (coeff(i) < 0)) )
								% -Var*coeff < 0  =>  z_i falls.
								% lower bound distance is important. 
		if z_b(i) == l(i)       % z_i would turn below l.
		  move = 0;
		  j_alternative = [j_alternative i];
		  xx_alternative = [xx_alternative i];
		  min_vec(i) = 0;
     	  %display:  
		  ds_zl=1;
		  Z_L_and_FALLS = [Z_L_and_FALLS i];%
		else
		  min_vec(i) = (z_b(i) - l(i)) / coeff(i); 
								% weighted distance to lower bound
		end

%%%%%%%%%%%%%%%%    upper bound    %%%%%%%%%%%%%%%%%%%%%%%%        

	  elseif (((direction == 1) && (coeff(i) < 0)) || ((direction == -1) && (0 < coeff(i))) )
								% -Var*coeff > 0  =>  z_i rises.
								% upper bound distance is important. 
		if z_b(i) == u(i)       % z_i would turn above u.  
		  move = 0;
		  j_alternative = [j_alternative i];
		  xx_alternative = [xx_alternative i];
		  min_vec(i) = 0;
     	  %display:  
		  ds_zu=1;
		  Z_U_and_RISES = [Z_U_and_RISES i];%		  
		else
		  min_vec(i) = (z_b(i) - u(i)) / coeff(i);  
								% weighted distance to upper bound
		end 
	  end  

%%%%%%%%%%%%%%%%    if coeff==0:    %%%%%%%%%%%%%%%%%%%%%%%%        

	else                        % z_i will not change in any case
	  min_vec(i) = ray;         
	end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  elseif (zwvt(i) == 2)         % ...each w variables in the basis
    if coeff(i) ~= 0

	  if (((direction == 1) & (0 < coeff(i))) | ((direction == -1) & (coeff(i) < 0)) )
								% -Var*coeff < 0  =>  w_i falls.
								% distance to zero is important. 
		if w_k(i) == 0			% w_i would turn below l.
		  move = 0;
		  j_alternative = [j_alternative i];
		  xx_alternative = [xx_alternative i];
		  min_vec(i) = 0;
       	  %display:  
		  ds_w=1;
		  W_ZERO_and_FALLS = [W_ZERO_and_FALLS i];
		else
		  min_vec(i) = w_k(i) / coeff(i); 
								% weighted distance to lower bound
		end
	  else                      % w_i rises
		min_vec(i) = ray;         
	  end

    %%%%%%%%%%%%%%%%    if coeff==0:    %%%%%%%%%%%%%%%%%%%%%%%%        
	else                        % w_i will not change in any case
	  min_vec(i) = ray;         
	  
	end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  elseif (zwvt(i) == 3)         % ...each v variables in the basis
    if coeff(i) ~= 0

	  if (((direction == 1) & (0 < coeff(i))) |  ((direction == -1) & (coeff(i) < 0)) )
								% -Var*coeff < 0  =>  v_i falls.
								% distance to zero is important. 
		if v_k(i) == 0			% v_i would turn below l.
					
		  move = 0;
		  j_alternative = [j_alternative i];
		  xx_alternative = [xx_alternative i];
		  min_vec(i) = 0;
     	  %display:  
		  ds_v=1;
		  V_ZERO_and_FALLS = [V_ZERO_and_FALLS i];%
		else
		  min_vec(i) = v_k(i) / coeff(i); 
								% weighted distance to lower bound
		end
	  else 
		min_vec(i) = ray;         
	  end

    %%%%%%%%%%%%%%%%    if coeff==0:    %%%%%%%%%%%%%%%%%%%%%%%%        
	else
 	  min_vec(i) = ray;          % v_i will not change in any case
	end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  elseif (zwvt(i) == 0)         % ...the t variable in the basis
								% should be in pivotupdate only!
	if minor == 0
	  if display_depth == 2
		disp('ERROR  : zwvt(i)=0 in pivotinit');
		i
	  end
	else
      if abs(coeff(i)) > 1e-12              %(i.e.   '~= 0'   )
		min_vec(i) = (t-1)/coeff(i);
		if ((direction == 1) & (min_vec(i) < 0))
		  min_vec(i) = ray;
		elseif ((direction == -1) & (min_vec(i) > 0 ))
		  min_vec(i) = ray;		
		end
	  else 
		min_vec(i) = ray;
	  end
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  elseif (zwvt(i) == 4)         % ...each a variables in the basis
    if coeff(i) ~= 0

	  if (((direction == 1) & (0 < coeff(i))) | ((direction == -1) & (coeff(i) < 0)) )
								% -Var*coeff < 0  =>  a_i falls.
								% distance to zero is important. 
		if a_k(i) == 0			% a_i would turn negative.
		  move = 0;
		  a_alternative = [a_alternative i];
		  xx_alternative = [xx_alternative i];
		  min_vec(i) = 0;% need this???
     	  %display:  
		  ds_a=1;
		  A_ZERO_and_FALLS = [A_ZERO_and_FALLS i];%
		else 
		  min_vec(i) = a_k(i) / coeff(i); 
								% weighted distance to zero
		end
	  else                      % a_i rises
		min_vec(i) = ray;         
	  end

    %%%%%%%%%%%%%%%%    if coeff==0:    %%%%%%%%%%%%%%%%%%%%%%%%        
	else                        % a_i will not change in any case
	  min_vec(i) = ray;         
	end 
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% In the following cases, the procedure does not move: 

% % display purposes:
% if display_depth == 2
%   if ds_v == 1
% 	V_ZERO_and_FALLS = V_ZERO_and_FALLS'
%   end
%   if ds_w == 1
% 	W_ZERO_and_FALLS = W_ZERO_and_FALLS'
%   end
%   if ds_a == 1
% 	A_ZERO_and_FALLS = A_ZERO_and_FALLS'
%   end
%   if ds_zu == 1
% 	Z_U_and_RISES = Z_U_and_RISES'
%   end
%   if ds_zl == 1
% 	Z_L_and_FALLS = Z_L_and_FALLS'
%   end
% end




% PARENTFILEs: pivotinit.m and pivotupdate.m