%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: vis3Dfield_cube
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
% Description: This function plots the flow field simulated in
% three-dimensions by the Stokes-Flow-Simulation. It additionally plots the
% boundaries elements of the surface, and shades the portion
% of the bottom surface that is moving. A 2D quiver plot is generated on
% the moving surface to indicate the magnitude and direction of flow
% 
% 
%
% Inputs:
% x, y, z: x and y coordinates of interior points
% u, v, w: flow vector of interior points
% bcx, bcy, bcz: x, y, z coordinates of boundary elements
% bcu, bcv, bcw, vbc: flow vector of boundary elements
% vertpts, elemverts: vertices and indices of triangles to plot of cube
% surface
% 
% Note, this script calls the plotsurface function from: 
% 
% as well as the quiver3D function from:
%
function vis3Dfield_cube(corner1,corner2,x,y,z,u,v,w,bcx,bcy,bcz,bcu,bcv,bcw,vertpts,elemvert)

% calculate some geometric properties of where the flow is calculated
xmin=corner1(1);
ymin=corner1(2);
zmin=corner1(3);
xmax=corner2(1);
ymax=corner2(2);
zmax=corner2(3);

cubevert=[xmin,ymin,zmin;xmin,ymin,zmax;xmin,ymax,zmin;xmin,ymax,zmax;...
          xmax,ymin,zmin;xmax,ymin,zmax;xmax,ymax,zmin;xmax,ymax,zmax];    
cubevert1=[cubevert(1,:);cubevert(2,:);cubevert(4,:);cubevert(3,:);cubevert(1,:); ...
           cubevert(5,:);cubevert(6,:);cubevert(8,:);cubevert(7,:);cubevert(5,:)];
cubevert2=[cubevert(1,:);cubevert(5,:);cubevert(6,:);cubevert(2,:)];
cubevert3=[cubevert(3,:);cubevert(7,:);cubevert(8,:);cubevert(4,:)];
         

bctot_ordered=sqrt(bcu.^2+bcv.^2+bcw.^2);
bcnonzero=bctot_ordered~=0;
bcxnz=bcx(bcnonzero);
bcynz=bcy(bcnonzero);
bcunz=bcu(bcnonzero);
bcvnz=bcv(bcnonzero);

% bcz and bcw here are not used because  we are doing a 2D quiver plot on
% the bottom surface.

% here we put the vertices and faces in a structure that is called by the
% plotsurface function
wallstl=struct('faces',elemvert,'vertices',vertpts);

% we can specify the color of the moving portion of the wall
coloring1=ones(size(elemvert,1),1)*0.5;
coloring1(bctot_ordered(:,1)~=0)=3;

% here we put all our positional information into a long array
posarray=[x(:), y(:), z(:)];

% here we calculate our speed so we can use the max to normalize our arrows
speedarray=sqrt(u(:).^2+v(:).^2+w(:).^2);

% if we just want to display the direction of the arrows without encoding
% total speed in the size, we can use normalize all speeds
% magarray=[u(:)./speedarray, v(:)./speedarray, w(:)./speedarray];

% Here, we set the magnitudes of the arrow size to be proprotional to the
% total speed
magarray=[u(:), v(:), w(:)]/3;

% Here we will also encode speed with color
dc=dopplerColors(64);
colarray=dc(floor(speedarray/max(speedarray(:))*32)+32,:);
n1=norm(corner1-corner2);

figure(30);
% call the plot surface function to plot the boundary elements
plotsurface(wallstl,coloring1,0.2);hold on;
% plot the edges of the cube for clarity
plot3(cubevert1(:,1),cubevert1(:,2),cubevert1(:,3),...
    cubevert2(:,1),cubevert2(:,2),cubevert2(:,3),...
    cubevert3(:,1),cubevert3(:,2),cubevert3(:,3),'Linewidth',3,'Color',[0,0,0])

% call the quiver3D function
quiver3D(posarray,magarray/max(speedarray(:))*n1,colarray);

% plot the 2D quiver field of the moving surface
quiver(bcxnz,bcynz,bcunz,bcvnz);axis off;
hold off;
camlight('headlight')

% Depending on your computer, you may have issues with this rendering. You
% can also try to not set this renderer, try 'zbuffer' (deprecated), or try
% 'painters'
set(gcf,'Renderer','opengl');
drawnow;