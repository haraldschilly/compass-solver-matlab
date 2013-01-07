%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPASS_driver: calls the COMPASS routine to solve a specified MCP.  
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

% This file illustrates the use of COMPASSmain.m by showing how to solve
% a mixed complementarity problem (MCP)  
 
% given  f: R^n -> R^n, and (possibly infinite) bounds l,u 
% find   z in [l,u], w,v in R_+^n
% s.t.:  f(z) = w-v
%        <(z-l),w> = 0
%        <(u-z),v> = 0

% The problem can be defined by the user, main options can be set, and
% the COMPASS solver is called: 

% The MCP is reformulated in terms of the normal equation, a pointwise
% general first order approximation (at the current iterate) of which is
% recast as a linear MCP and tried to be solved in a pivot procedure.
% The point returned by this pivot procedure is taken to be the next
% iterate.  
% A non-monototne stabilization scheme in combination with a watchdog
% technique is used in the outer algorithm so as to minimize a smooth
% merit function, the zero of which is known to be equal to a solution
% of the MCP. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create problem/ problem definition:

n=10;                    % problem dimension       
x0= unifrnd(-1,2,1,n);   % initial value
l = unifrnd(-1,1,1,n);   % upper and lower bounds 
u = l + rand(1,n);       

% There are two ways to create the objective function f(x) and Jacobian
% Df(x) (which is a necessary input from the user): 

% 1. The expression of the objective function can be stated in the file
% objfcn.m, that of the Jacobian in the file Jacobian_objfcn.m. 

% f = @objfcn;
% Df = @Jacobian_objfcn;      

% 2. The objectivefunction and Jacobian can be defined here directly, by
% typing the expressions after f=@(x) and Df=@(x), respectively.  

% f=@(x)   <your expression here> ;
% Df=@(x)  <your expression here> ;

% To solve your own problem, simply change the bounds and function
% expressions above. Be sure to have '%' signs in front of lines that
% should NOT be read by octave/MATLAB.  


% Example: a randomly generated general quadratic function: 
 DM= diag(unifrnd(-1,1,1,n));
 A = unifrnd(-1,1,n);
 b = unifrnd(-1,1,1,n);
 f=@(x) (DM*x'.^2 + A*x' + b')';
 Df=@(x) 2.*diag(DM*x') + A;




%%%%%%%%%%%%%%%%%%%%%
% main options: 

max_major=15*n;        % restart the algorithm after max_majoriter
                       % major iterations.  
max_minor=15*n;        % stop the pivot procedure after
                       % max_minoriter minor iterations. 
max_restarts=20;       % stop COMPASS after min12 to max60 restarts 
                       % without a solution. 
target=1.e-8;          % quit when a point with merit function value
                       % <= target is found
display_depth = 0;     % display details on solving process: 
					   % 0 -solution only 
					   % 1 -major iteration info 
					   % 2 -major and minor iteration info



%%%%%%%%%%%%%%%%%%%%%%
% call solver:
disp(' ');
disp('         solving by COMPASS ...');disp(' ');
tic;

COMPASSmain

disp(' ');
disp('          ... solved by COMPASS:');disp(' ');
duration_in_seconds = toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display output:
if Psi_value < target
  disp(' ');
  disp('     l         z_k         u          w_k        v_k')
  Solution = [l' z_k' u' w_k' v_k']  
  disp(' ');
  disp('call, major_iterations, av.minor_iterations, last_Psi')
  display_statistics = [(1:restarts+1)' tot_major' save_avg_minor' save_Psi']
else
  disp('no solution found');
  disp(' ');
  disp('     l         z_k         u          w_k        v_k')
  last_points = [l' z_k' u' w_k' v_k']  
  disp(' ');
  disp('display_statistics:');
  disp('call, major_iterations, av.minor_iterations, last_Psi') 
  [(1:restarts+1)' tot_major' save_avg_minor' save_Psi']
end
Psi_value
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ');
disp(' ');



