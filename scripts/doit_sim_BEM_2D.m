%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: doit_sim_BEM_2D
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
% cavity using the boundary element method (BEM).
%
% 
% Inputs: 
% BCType - type of boundary conditions, either traction or velocity
% boundundbound - set as bound or unbound. The unbound (open) geometry consists
% of only a single horizontal surface, with fluid free to circulate
% elsewhere. The bound (closed) geometry is a rectangular cavity
% corner1, corner2 - these two points need to be 2d coordinates that give
% the lower left and upper right boundary of the rectangular cavity. In the
% case of the unbound condition, they will also define where we calculate
% our flow
% nbc - number of boundary elements along each of the four walls of the
% surface (total points will be nbc x 4 for closed geometry)
% Ubound - set the speed (or traction) of the moving portion of the boundary
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
% ub,vb: flow vector of boundary elements
% fx,fy: traction at boundary elements
%
%  
%% Specify the inputs as detailed in the script description
corner1=[0,0];
corner2=[0.25,0.125];
xpts=20;
ypts=10;
numbc=40;
Uspeed=10;
fos=0.8;
visc=1;
%% Generate the geometry and boundary conditions of the flow
[elemvert,vertpts,bcx,bcy,bcu,bcv,p,normv,alpha] = generateBC_BEM_2D('traction','bounded',corner1,corner2,numbc,Uspeed,fos);
%% Plot the location of the boundary points as well as the

dx=corner2(1)-corner1(1);
dy=corner2(2)-corner1(2);
figure(1);plot(p(:,1),p(:,2),'x',vertpts(:,1),vertpts(:,2))
xlim([corner1(1)-0.1*dx corner2(1)+0.1*dy]);
ylim([corner1(2)-0.1*dy corner2(2)+0.1*dy]);
title('Plot of boundaries of problems, with x marking center of each boundary element')
xlabel('x coord')
ylabel('y coord')
%% Calculate the Green's Function matrix (single and double layer potential)
[SLmat,DLmat]=calc_Amat_BEM_2D(p,vertpts,elemvert,normv,true);
%% Solve the boundary integral equation (see tutorial Eq. 41)
BCtotal=[bcu;bcv];
[u0,f0,ub,vb,fx,fy] = solve_mixedBC(SLmat,DLmat,BCtotal,visc,alpha,'2D');

%% Generate a grid over the points we would like to evaluate in the domain
% High shearing occurs near the corners, so numerical errors are increased
% in those locations. The last parameter gives how close to these corners
% we would like to calculate
[xint,yint]=InteriorPts_2D(corner1,corner2,xpts,ypts,0.95);
intpts=[xint(:),yint(:)];

%% Calculate the single and double layer potential effects on the interior points
[SLint,DLint]=calc_Amat_BEM_2D(intpts,vertpts,elemvert,normv,false);

%% Calculate the pressure and shear stress single and double layer potentials
[PSLint,PDLint,SSLint,SDLint]=calc_PSmat_BEM_2D(intpts,vertpts,elemvert,normv,false);

%% Calculate the flow field u, v, option to also calculate the pressure and shear stress tensors
%[u,v] = calcUV_BEM(f0,u0,SLint,DLint,xpts,ypts,visc);
%[u,v,pres] = calcUV_BEM(f0,u0,SLint,DLint,xpts,ypts,visc,PSLint,PDLint);
[u,v,pres,stress]=calcUV_BEM_2D(f0,u0,SLint,DLint,xpts,ypts,visc,PSLint,PDLint,SSLint,SDLint);

%% Plot the flow field and boundary conditions
vis2Dfield_quiver(xint,yint,u,v,bcx,bcy,ub,vb,2);