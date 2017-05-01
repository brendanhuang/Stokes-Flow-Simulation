%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: InteriorPts_3D
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
% Description: This function generates a grid of points on the interior of
% the domain of our Stokes flow problem in two dimensions.
% 
%
% Inputs:
% corner1, corner2 - corners of cuboid to define the meshgrid
% xpts, ypts, zpts - number of points in the x, y, and z dimension
% closeness - how close we want to put points to the corners, where we may
% encounter high shearing
%
% Outputs: 
% xc, yc, zc - cubic grid of x, y, and z coordinates

function [xc,yc,zc]=InteriorPts_3D(corner1,corner2,xpts,ypts,zpts,closeness)

dx=(corner2(1)-corner1(1))*(1-closeness);
dy=(corner2(2)-corner1(2))*(1-closeness);
dz=(corner2(3)-corner1(3))*(1-closeness);

[xc,yc,zc]=ndgrid(linspace(corner1(1)+dx,corner2(1)-dx,xpts),linspace(corner1(2)+dy,corner2(2)-dy,ypts),linspace(corner1(3)+dz,corner2(3)-dz,zpts));
end