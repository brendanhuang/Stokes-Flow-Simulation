%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: generateBC_BEM_2D
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
% Description: This function generates the boundary points and boundary
% conditions for the cavity driven fluid flow. There are several options
% for input
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
%
% Outputs:
% vertpts - location of vertices of line (2D) of boundaries 
% elemvert - vector to identify which vertices belong to which boundary element
% p - location of points where boundary elements are evaluated;
% bcx, bcy, bcu, bcv - values of boundary element points, both location and value of either velocity or traction; 
% alpha - index specifying at each of bcu, bcv, whether it is a velocity or traction condition
% normv - normal vector of each boundary element

function [elemvert,vertpts,bcmidx,bcmidy,bcu,bcv,p,normv,alpha] = generateBC_BEM_2D(BCtype,boundunbound,corner1,corner2,nbc,Ubound,fos)

% in case we don't specify the fraction of the moving surface, we'll assume
% the whole thing is moving
if nargin==6
    fos = 1;
end

% calculate the total width and height of the cavity
dx=corner2(1)-corner1(1);
dy=corner2(2)-corner1(2);

% generate the location of the boundary elements and the speed / traction
% associated with each
% in the case of unbounded, we only need to worry about the bottom surface
% in the case of bounded cavity, we need to construct four walls
switch boundunbound
    case 'unbounded'
        bxbot = linspace(corner1(1),corner2(1),nbc)';
        bybot = linspace(corner1(2),corner1(2),nbc)';
        bcu= linspace(0,0,nbc)';
        bcv= linspace(0,0,nbc)';
        velboundleft=round(nbc/2-fos*nbc/2);
        velboundright=min(round(nbc/2+fos*nbc/2),nbc-1);
        bcu(velboundleft:velboundright)=Ubound;
        bcxcoor=(bxbot(1:end-1)+bxbot(2:end))/2;
        bcycoor=(bybot(1:end-1)+bybot(2:end))/2;
        bcu=bcu(1:end-1);
        bcv=bcv(1:end-1);
        l1=linspace(1,nbc,nbc)';
        elemvert=[l1,circshift(l1,[1,1])];
        elemvert=elemvert(2:end,:);
        vertpts=[bxbot,bybot];
    case 'bounded'
        bxbot=linspace(corner1(1),corner2(1),nbc)';
        bybot=linspace(corner1(2),corner1(2),nbc)';
        bxright=linspace(corner2(1),corner2(1),nbc)';
        byright=linspace(corner1(2),corner2(2),nbc)';
        bxtop=linspace(corner2(1),corner1(1),nbc)';
        bytop=linspace(corner2(2),corner2(2),nbc)';
        bxleft=linspace(corner1(1),corner1(1),nbc)';
        byleft=linspace(corner2(2),corner1(2),nbc)';
        bxbot=bxbot(1:end-1);bybot=bybot(1:end-1);
        bxright=bxright(1:end-1);byright=byright(1:end-1);
        bxtop=bxtop(1:end-1);bytop=bytop(1:end-1);
        bxleft=bxleft(1:end-1);byleft=byleft(1:end-1);
        bcxcoor=[bxbot; bxright; bxtop; bxleft];
        bcycoor=[bybot; byright; bytop; byleft];        
        
        bcv=zeros(size(bcxcoor));
        velboundleft=max(round(nbc*(1/2-fos/2)),2);
        velboundright=min(round(nbc*(1/2+fos/2)),nbc);
        bubot=linspace(0,0,nbc)';
        bubot(velboundleft:velboundright)=Ubound;
        buright=linspace(0,0,nbc)';butop=linspace(0,0,nbc)';buleft=linspace(0,0,nbc)';
        bubot=bubot(1:end-1);buright=buright(1:end-1);butop=butop(1:end-1);buleft=buleft(1:end-1);
        bcu=[bubot; buright; butop; buleft];
        l1=linspace(1,numel(bcxcoor),numel(bcxcoor))';
        elemvert=[l1,circshift(l1,[1,1])];
        vertpts=[bcxcoor,bcycoor];       
end

% qa and qb are intermediate variables to help us calculate p, the 
% midpoints of each boundary element. Note that when we calculate our
% initial matrix, we calculate the effect of a full element (line segment)
% on a single point that is taken to be the midpoint of another line
% element
qa=[vertpts(elemvert(:,1),1), vertpts(elemvert(:,1),2)];
qb=[vertpts(elemvert(:,2),1), vertpts(elemvert(:,2),2)];

p=(qa+qb)/2;
bcmidx=p(:,1);
bcmidy=p(:,2);
   
% calculate a unit vector normal to each element.
pm=qa-qb;
qlen=sqrt(pm(:,1).^2+pm(:,2).^2);
normv=zeros(size(pm));
normv(:,1)=-pm(:,2)./qlen; 
normv(:,2)=pm(:,1)./qlen; 

% now we calculate alpha, which assigns a value to each boundary element 
% stating whether that element specifies a velocity or traction boundary
% condition. Notably, we may specify the y-direction as flow (impenetrable
% wall), but specify the x-direction as traction
switch BCtype
    case 'velocity'
        alpha = false(size(bcu));
        alpha = [alpha;alpha];
    case 'traction'
        alpha = zeros(size(bcu));
        alpha(bcu~=0)=1;
        alpha=logical([alpha;zeros(size(bcu))]);
        bcu=-bcu;
end
    
end