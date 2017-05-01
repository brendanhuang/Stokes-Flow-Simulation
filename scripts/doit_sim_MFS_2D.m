%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: doit_sim_MFS_2D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MIT License
% 
% Copyright (c) 2015 Brendan K. Huang
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description: This doit script runs a simulation of a 2D lid-driven
% cavity using the method of fundamental solutions (MFS).
%
% 
% Inputs: 
% boundundbound - set as bounded or unbounded. The unbound (open) geometry consists
% of only a single horizontal surface, with fluid free to circulate
% elsewhere. The bound (closed) geometry is a rectangular cavity
% corner1, corner2 - these two points need to be 2d coordinates that give
% the lower left and upper right boundary of the rectangular cavity. In the
% case of the unbound condition, they will also define where we calculate
% our flow
% nbc - number of boundary points / stokeslets along each of the four walls of the
% surface (total points will be nbc x 4 for closed geometry) that we will
% evalute
% Ubound - set the speed of the moving portion of the boundary. Traction
% currently not implemented
% fos - fraction of surface of bottom wall that is moving or applying
% traction to the fluid. fos = 1 means the entire wall moves
% xpts, ypts: number of points in the x and y direction to be evaluated on
% the interior of the domain for flow
% visc - the viscosity of the fluid is required to relate traction to flow
%
% Outputs:
% xint,yint: grid of x,y coordinates where flow is evaluated
% u,v: flow vector at locations xint, yint
% bcx,bcy: location of centers of boundary elements
% bcu,bcv: flow vector of boundary elements
%
%  
%% Specify the inputs as detailed in the script description
corner1=[0,0];
corner2=[1,1];
xpts=20;
ypts=20;
numbc=20;
b=0.01;
Uspeed=10;
fos=0.5;
%% Generate the geometry and boundary conditions of the flow
[scx,scy,bcx,bcy,bcu,bcv]=generateBC_MFS_2D(1,'bounded',corner1,corner2,numbc,b,Uspeed,fos);

%% Calculate the Green's Function matrix 
Amat=calc_Amat_MFS_2D(bcx,bcy,scx,scy);
BCtotal=transpose([bcu,bcv]);

%% Solve the inverse equation (see tutorial Eq. 14)
qcoeff=inv(Amat)*BCtotal;
%% Generate a grid over the points we would like to evaluate in the domain
% High shearing occurs near the corners, so numerical errors are increased
% in those locations. The last parameter gives how close to these corners
% we would like to calculate
[xint,yint]=InteriorPts_2D(corner1,corner2,xpts,ypts,0.95);

%% Calculate the flow, pressure, and stream function Green's functions on the interior
[Aint,Pint,Psiint]=calc_Amat_MFS_2D(transpose(xint(:)),transpose(yint(:)),scx,scy);

%% Compute the flow field, pressure, and stream function
[u,v] = calcUV_MFS_2D(qcoeff,Aint,xpts,ypts);
press = calcPPsi_MFS_2D(qcoeff,Pint,xpts,ypts);
psi= calcPPsi_MFS_2D(qcoeff,Psiint,xpts,ypts);

%% Plot the flow field, streamlines, and boundary conditions
vis2Dfield_quiver(xint,yint,u,v,bcx,bcy,bcu,bcv,1);
vis2Dfield_streamline(xint,yint,psi,bcx,bcy,bcu,bcv,2)
