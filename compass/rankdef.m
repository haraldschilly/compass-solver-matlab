%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rankdef: Curing rank deficiency by inserting artificial variables. 
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

% Replace lin.dep. 'z'-kind columns in Basis by artificial 'a' kind
% variables and columns (unity matrix columns of the respective indices)
% so that Basis is invertible. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Basis_new = Basis_old;
H_new = H_old;
b_vec_new = b_vec_old;

d_ = zeros(n);
delta = zeros(1,n);

% Replace the z-kind column that greatest
II=eye(n);
for i = 1:n
  if zwvt(i)==1
	d_(:,i) = II(:,i)-Basis(:,i);	% difference vectors (replace one!)  
	delta(i) = (1+H_old(i,:)*d_(:,i));  % denominator of each possible z-variable
  end
end

[dummy, ldc] = max(abs(delta)); % replaceable index: ldc (linear dependent column)   

if dummy > 1e-8
  delta_ldc = delta(ldc);       % the value of delta for the replaceable index...

  % replace:
  Ia(ldc,ldc) = 1;
  M_rem(:,ldc) = M(:,ldc);
  z_rem(ldc) = z_b(ldc);        % now in M_rem but not at bound! 
  M_b(:,ldc) = zeros(1,n)';                  
  z_b(ldc) = 0;
  zwvt(ldc) = 4;                % a kind columns have zwvt-index 4. 

  % update inverse (H), Basis and b_vec:
  beta = 1/delta_ldc;
  H_new = H_old - beta.* H_old*d_(:,ldc) * H_old(ldc,:);
  b_vec_new = z_b + w_k + v_k + a_k;
  Basis_new = M_b - Iw + Iv + Ia;   

  if abs(det(Basis_new*H_new)-1) > 1e-8
	if display_depth == 2
	  disp(' det(Basis_new*H_new) ~= 1: abort ');
	end
	abort = 1;
	unsuccessful = 1;
 	take_pg_step = 1;	 
  end
  
else
  if display_depth == 2
	disp('Basis rank deficient. No cure possible. ');
  end
  abort = 1;
  unsuccessful = 1;
  take_pg_step = 1;
%	this can have 3 reasons: 
%	1. r causes the problem: there is no z in the basis, the replacing of which
%	would cure the problem.
%	2. No replacing of any z would yield a new STABLE invertible Basis.
%   3. THere is no z in the Basis. -> r causes the problem.	(d_ = delta = 0)
end



% PARENTFILES: pivotinit.m and pivotupdate.m