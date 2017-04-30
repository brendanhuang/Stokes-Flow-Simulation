%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: vis2Dfield_quiver
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
% two-dimensions by the Stokes-Flow-Simulation. It additionally plots the
% boundary points and the flow on the boundary. Square represent stationary
% points on the boundary
% 
% 
%
% Inputs:
% x, y: x and y coordinates of interior points
% u, v: flow vector of interior points
% xbc, ybc: x and y coordinates of boundary elements
% ubc, vbc: flow vector of boundary elements
% fignum: optional to specify numbering of Matlab figure
% 

function vis2Dfield_quiver(x,y,u,v,xbc,ybc,ubc,vbc,fignum)

if nargin<9
    fignum=1;
end

% Find elements where there is flow on the boundary
ind1=bitor(ubc~=0,vbc~=0);
if numel(ubc(ind1))>10
% for visualization sake, we don't plot the flow at every boundary point if
% we have more than 10. We'll use our function reducemat to downsample our
% points
   % ds is downsample factor
   ds=floor(numel(ubc(ind1))/10); 
   
   % leftover is how many values we omit at the end
   leftover=mod(numel(ubc(ind1)),ds);
   unz=reducemat(ubc(ind1),ds,leftover);
   vnz=reducemat(vbc(ind1),ds,leftover);
   xnz=reducemat(xbc(ind1),ds,leftover);
   ynz=reducemat(ybc(ind1),ds,leftover);
else
   unz=ubc(ind1);
   vnz=vbc(ind1);
   xnz=xbc(ind1);
   ynz=ybc(ind1);
end

% Find the boundaries of domain
x1=min(x(:));
x2=max(x(:));
y1=min(y(:));
y2=max(y(:));
xl=x2-x1;
yl=y2-y1;

% set color of boundary points
boundarycolor=[0 0 0]/255; 
figure(fignum);

% call the color quiver function
quiverc(x(:),y(:),u(:),v(:));hold on;

% plot the boundaries
plot(xbc(~ind1),ybc(~ind1),'s','MarkerEdgeColor',boundarycolor,'LineWidth',1.5);
axis([x1-0.1*xl, x2+0.1*xl, y1-0.1*yl, y2+0.1*yl]);

% quiver the moving portion of the boundary, but rescale by physical size
% of area and max speed of moving wall
speed=sqrt(unz.^2+vnz.^2);
fact=max(speed)/xl;
g1=quiver(xnz,ynz,unz,vnz,0.3/fact,'Linewidth',1.5);
set(g1, 'Color',boundarycolor);
hold off;
end

% reducemat a function to downsample the quiver field if there are too many
% points on the boundary
function zz=reducemat(mat1,ds,leftover)
   a2=mat1(1:end-leftover);
   a3=mat1(end-leftover:end);
   a2=squeeze(mean(reshape(a2,ds,[]),1));
   a3=mean(a3);
   zz=[a2,a3];
end