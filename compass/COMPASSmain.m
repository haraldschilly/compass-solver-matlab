%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPASSmain: Main file of the COMPASS routine
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

% COMPASS is an algorithm to solve the mixed complementarity problem
% (MCP):

% given  f: R^n -> R^n, and (possibly infinite) bounds l,u 
% find   z in [l,u], w,v in R_+^n
% s.t.:  f(z) = w-v
%        <(z-l),w> = 0
%        <(u-z),v> = 0

% A crash technique is used at the beginning of the algorithm in order
% to determine an approximation to the active set at the solution. 

% Solving the normal equation, a reformulation of the MCP in terms of
% the nonsmooth normal map, as proposed by Robinson, is the core of the
% algorithm. A Newton like method is used:

% A general first order approximation to the normal map is
% computed. The problem of finding a zero of this approximation is
% recast as a linear MCP. This linear MCP is solved by a pivot
% technique similar to that of Lemke, or that described by Dantzig.  
% The pivot procedure (inner, or minor algorithm) yields a piecewise
% linear path connecting the current iterate and the Newton point. Each
% pivot step results in a linear piece of the path. 

% If the new iterate does not provide the necessary descent in the
% merit function, the constructed path (or the line segment connecting
% its endpoints) is searched for an appropriate point.  
% A Non monotone stabilization scheme (NMS) and a watchdog techniqe are
% used in order to reduce the number of function and gradient
% evaluations necessary (outer, or major algorithm).  

% This procedure is repeated until the normal equation, and hence the
% MCP, is solved. 

% An example for how to use the program is in COMPASS_driver.m

% The program was written by Stefan Schmelzer (University of Vienna) as
% part of his diploma thesis. The current Version 1.0 is dated from
% April 16, 2012. 
 
% Please inform the author at schmelzer@ihs.ac.at if you make 
% serious use of this code. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%
% Options:            % (can be changed by user; some are automatically
                      % changed after a restart)   

% crash technique
do_crash = 1;           % flag for running/stopping the crash procedure 
min_dim_crash = 10;     % minimum problem dimension: n=10
sigma_crash = 0.05;     % minimum descent in crash technique  
max_m_alpha_crash = 14; % minimum step size in crash technique is
                        % (0.5)^12   
max_crash_count = 50;   % maximum number of crash iterateions
min_active_change = 1;  % minimum change in active index set 
max_no_change_row = 3;  % max number of consecutive steps with
						% insufficient change in active set.  

% pivot technique
in_minoriter = 0;      % flag: is the algorithm in the file minoriter? 
in_pivot = 0;		   % flag: is the algorithm in the pivot procedure?  
infty=1000;            % number used as infinity border
cyclesize = min(max(5, n/10), 20);   % size that the length of a cycle
               						 % (of pivotroutine) is assumed.   
maybe_cycle = 2;       % iteration at which cycling possibility is
					   % checked.   
definitive_cycling = 0.5*n;   % maximum number of the same variable
							  % entering without abort due to cycling   
max_no_solution_LCP_row = 10; % maximum number of unsuccessful pivot
							  % procedures in a row before a restart.   

% NMS
pathsearch_type = 1;   % chose: line(1) or arc(2)
pg_search_type = 2;    % chose: line(1) or arc(2)
sigma_NMS = 0.01;      % minimum descent in pg and watchdog searches
max_m_alpha_w = 20;    % minimum step size in watchdog searches
max_m_alpha_pg = 25;   % minimum step size in pg searches
max_d_count = 3;       % maximum number of consecutive d-steps 
Delta0 = 0.5;          % measure of closeness in every first of
                       % consecutive d-steps  
deltachange = 0.07;    % Delta is reduced by deltachange after
                       % consecutive d-steps 

% restarts: 
next_x_near = 0;       % flag to take a special x0 after restart 
x_next = norm(u-l)*rand(1,n); % default special x0 after restarts 
max_restarts = min(60,max(12,max_restarts));
					   % set max_restarts to min.12 and max.60 restarts. 

% statistics:
restarts = 0;          % number of restarts in the algorithm
tot_major = [];        % major iterations per model run (before a 
                       % restart is made)  
save_avg_minor = [];   % average pivot steps per model run
save_Psi =[];          % last Psi_value before restart or total stop.




% initial settings: 					
x_k = x0;                            % current point
z_k = piB(x_k,l,u);                  % projection onto box B 
Psi_value_0 = Psi(z_k,l,u,f,infty);  % merit fct value at x0.
Psi_value_l = Psi(l,l,u,f,infty);    % merit fct value at l.
Psi_value_u = Psi(u,l,u,f,infty);    % merit fct value at u.
Psi_value_init = Psi_value_0;
if display_depth ~= 0
  Psi_value_0
  Psi_value_init
end


%%%%%%%%%%%%%%%%%%
% call the crash technique:
crash            
% Computes an approximation to the active set at the solution. 
% Returns z_k, x_k, w_k and v_k
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%
% call the core of the algorithm
majoriter
% iteratively solves MCP
%%%%%%%%%%%%%%%%%%%



if restart == 1    
 % The algorithm could not solve the problem with the initial settings.
 % It is restarted with different options, multiple times if necessray.

  while restart == 1               

	%%%%%%%%%%%%%%%%
    % call the restart routine:
	RestartCOMPASS
	% different options are defined
	% crash technique is called
	% majoriter is called. 
	%%%%%%%%%%%%%%%%%%%%%

  end
end
