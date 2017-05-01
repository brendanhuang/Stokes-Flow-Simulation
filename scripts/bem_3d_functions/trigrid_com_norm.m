%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: trigrid_com_norm
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
% Description: This function takes the coordinates of a mesh grid and the
% output of the grid2tri function (specifying indices) and returns the
% center of mass as well as the normal vector of the triangles generated on
% the grid
% 
% Inputs: 
% tri - indexing output of tri2grid
% xcoor, ycoor, zcoor - meshgrid at all locations, split up into x, y, and
% z coordinates
%
% Outputs:
% normd - unit normal vector of each triangle
% com - center of mass of each triangle. The center of mass will be used as
% our points to evaluate on the each boundary element

function [normd,com] = trigrid_com_norm(tri,xcoor,ycoor,zcoor)
com=zeros(size(tri));

com(:,1)=mean([xcoor(tri(:,1)) , xcoor(tri(:,2)), xcoor(tri(:,3))],2);
com(:,2)=mean([ycoor(tri(:,1)) , ycoor(tri(:,2)), ycoor(tri(:,3))],2);
com(:,3)=mean([zcoor(tri(:,1)) , zcoor(tri(:,2)), zcoor(tri(:,3))],2);

qa=[xcoor(tri(:,1)),ycoor(tri(:,1)),zcoor(tri(:,1))];
qb=[xcoor(tri(:,2)),ycoor(tri(:,2)),zcoor(tri(:,2))];
qc=[xcoor(tri(:,3)),ycoor(tri(:,3)),zcoor(tri(:,3))];

u=qb-qa; 
v=qc-qa;
Nx=u(:,2).*v(:,3)-u(:,3).*v(:,2);
Ny=u(:,3).*v(:,1)-u(:,1).*v(:,3);
Nz=u(:,1).*v(:,2)-u(:,2).*v(:,1);
mag=sqrt(Nx.^2+Ny.^2+Nz.^2);
normd=[Nx./mag,Ny./mag,Nz./mag];



end