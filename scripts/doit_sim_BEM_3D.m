%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: doit_sim_BEM_3D
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
% Description: This doit script runs a simulation of a 3D lid-driven
% cavity using the boundary element method (BEM).
%
% 
% Inputs: 
% BCType - type of boundary conditions, either traction or velocity
% boundundbound - set as bound or unbound. The unbound (open) geometry consists
% of only a single horizontal surface, with fluid free to circulate
% elsewhere. The bound (closed) geometry is a rectangular cavity
% corner1, corner2 - these two points need to be 3D coordinates that give
% the bottom front left and top back right boundary of the cuboid cavity. In the
% case of the unbound condition, they will also define the limeits we calculate our flow
% nbc - number of boundary elements along each of the four walls of the
% surface (total points will be nbc x 6 for closed geometry)
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
%% Define basic parameters of simulation
corner1=[-5,-5,0]; % bottom left corner of cube where flow will be simulated
corner2=[5,5,5]; % top right corner of cube
xpts=6; % number of x, y, and z pts to simulate in interior
ypts=6;
zpts=3;
numbc=5; % size of mesh
Uspeed=1; % either speed or traction at boundary
fos=0.2; % fraction of surface that is "moving"
visc=1; % viscosity of the fluid

%% generate boundary conditions and geometry
% triangular mesh, consisting of vertpts (locations of triangle vertices)
% and elemvert (describes the vertices of each triangle), as well as
% boundary conditions (specifying 'velocity' or 'traction' [mixed]), center
% points of the triangle p, normal vectors at those locations, and the
% alpha vector, which denotes which boundary points are velocity, and which
% are traction
%
[elemvert,vertpts,bcx,bcy,bcz,bcu,bcv,bcw,p,normv,alpha] = generateBC_BEM_3D('traction','bounded',corner1,corner2,numbc,Uspeed,fos);
%% plot the mesh as well as the normal vectors to double check that everything looks kosher
% For the closed domain, the normal vectors should all be pointing inward
% to the cube.
figure(1);trimesh(elemvert,vertpts(:,1),vertpts(:,2),vertpts(:,3));
title('Mesh of boundary')
figure(2);quiver3(p(:,1),p(:,2),p(:,3),normv(:,1),normv(:,2),normv(:,3));
title('Normal vectors of each element, plotted at location of element')
%% calculate the single and double layer potential (i.e. susceptibility matrix)
% between each element on the boundary boundary.  Special algorithms are
% required when calculating susceptibility between an element and itself (singular),
% and ponq denotes that we are considering singular points
tic
[SLmat,DLmat]=calc_Amat_BEM_3D(p,vertpts,elemvert,normv,true);
toc
%% solve the boundary integral equation (see tutorial Eq. 41)
% solve_mixedBC will solve the Fredholm integral of the first kind for
% velocity. If the boundary conditions are mixed, it will swap out the
% unknown velocities for known stresses and will end up still solving the
% Fredholm integral.  Robin conditions are possible to solve but
% not addressed here. f0 is the force at the boundary (3 components each),
% while u0 is the velocity at the boundary (3 components each).
BCtotal=[bcu;bcv;bcw];
[u0,f0,ub,vb,wb,fx,fy,fz] = solve_mixedBC(SLmat,DLmat,BCtotal,visc,alpha,'3D');
%% define a 3D grid where you want to evaluate the fluid in the box

[xint,yint,zint]=InteriorPts_3D(corner1,corner2,xpts,ypts,zpts,0.90);
intpts=[xint(:),yint(:),zint(:)];
%% Calulate pressure and stress susceptibility matrices
% both single layer and double layer potentials for each
% PSL = pressure single layer, PDL = pressure double layer
% SSL = stress single layer, SDL = stress double layer
[PSLint,PDLint,SSLint,SDLint]=calc_PSmat_BEM_3D(intpts,vertpts,elemvert,normv,false);

%% Calulate single and double layer potential for flow for interior points
tic
[SLint,DLint]=calc_Amat_BEM_3D(intpts,vertpts,elemvert,normv,false);
toc
%% Plug these back into the boundary integral equations to evaluate
% flow field, pressure, and shear stress tensor at interior points
% because stress tensor is 3x3, it has been collapsed into a single tensor,
% while the 1x3 velocity vector is given as separate components

[u,v,w] = calcUV_BEM_3D(f0,u0,SLint,DLint,xpts,ypts,zpts,visc);
%[u,v,w,pres] = calcUV_BEM_3D(f0,u0,SLint,DLint,xpts,ypts,zpts,visc,PSLint,PDLint);
%[u,v,w,pres,stress] = calcUV_BEM_3D(f0,u0,SLint,DLint,xpts,ypts,zpts,visc,PSLint,PDLint,SSLint,SDLint);

%% Visualize flow field
vis3Dfield_cube(corner1,corner2,xint,yint,zint,u,v,w,bcx,bcy,bcz,ub,vb,wb,vertpts,elemvert)
