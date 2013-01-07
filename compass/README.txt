%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
readme: COMPASS
Copyright (C) 2012 Stefan Schmelzer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This is the readme file of COMPASS.

COMPASS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

COMPASS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
  
You should have received a copy of the GNU General Public License
along with COMPASS.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


COMPASS is an algorithm to solve the mixed complementarity problem
(MCP):

given  f: R^n -> R^n, and (possibly infinite) bounds l,u 
find   z in [l,u], w,v in R_+^n
s.t.:  f(z) = w-v
       <(z-l),w> = 0
       <(u-z),v> = 0

A crash technique is used at the beginning of the algorithm in order
to determine an approximation to the active set at the solution. 

Solving the normal equation, a reformulation of the MCP in terms of
the nonsmooth normal map, as proposed by Robinson, is the core of 
the algorithm. Finding a zero of the normal map is equivalent to 
solving the MCP. A Newton like method is used:   

A general first order approximation to the normal map is
computed. The problem of finding a zero of this approximation is
recast as a linear MCP. This linear MCP is solved by a pivot
technique similar to that of Lemke, or that described by Dantzig.  
The pivot procedure (inner, or minor algorithm) yields a unique 
point, called Newton point, which is the first candidate to be the 
next iterate, an a piecewise linear path connecting the current 
iterate and the Newton point. Each pivot step results in a linear 
piece of this path. 

For the outer algorithm, a smooth merit function is used. A zero of 
this merit function also corresponds to a solution of the MCP. 
If the Newton point does not provide the necessary descent in the
merit function, the constructed path (or the line segment connecting 
its endpoints) is searched for an appropriate point. However, 
descent in the merit function is not necessary in every step. The 
algorithm is allowed to accept points with an increase in the merit 
function. This is done with the help of a Non monotone stabilization 
scheme (NMS) and a watchdog techniqe, in order to reduce the number 
of function and gradient evaluations necessary (in the outer, or 
major algorithm).  

This procedure is repeated until the value of the merit function is 
lower than a specified number ("target"=1e-8). In this case a 
solution to the normal equation is found, and hence the MCP is 
solved. 

The algorithm has a strong global convergence theory, and typically 
yields descent during progress in the outer algorithm. However some 
heuristics are used in the form of restarts, e.g. when progress 
becomes too slow near a local minimizer of Psi. 

An example of how to use the program is given in COMPASS_driver.m. 
The starting point x0, the borders l and u, the function f as well as 
the Jacobian of the function f are necessary inputs the user has to 
provide in order to be able to solve her problem. 
A considerable number of additional settings can also be changed by 
the user in COMPASSmain.m in order to reach best performace for her 
specific problem. Performance depends on the combination of the used 
options and problem characteristics. 

The program was written by Stefan Schmelzer (University of Vienna) as
part of his diploma thesis. The current version 1.0 is dated from
April 16, 2012. It is compatible with MATLAB and GNU Octave. 
 
Please inform the author at schmelzer@ihs.ac.at if you make serious 
use of this code. 

