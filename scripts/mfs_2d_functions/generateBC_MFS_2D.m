%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: generateBC_MFS_2D
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
% conditions for the cavity driven fluid flow. 
% 
% Inputs: 
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
% b - how far away we would like to place the Stokeslets. In this script,
% stokeslets are placed along the line between a boundary point and the
% center of mass (com), at a distance dtermined by b
% fos - fraction of surface of bottom wall that is moving or applying
% traction to the fluid. fos = 1 means the entire wall moves
% ploton - whether to plot the boundary points to review
%
% Outputs:
% stokescoordx, stokesccordy - location of stokeslets (fictious forces)
% bcxcoor, bcycoor, bcu, bcv - values of boundary element points, both location and value of either velocity or traction; 

function [stokescoordx,stokescoordy,bcxcoor,bcycoor,bcu,bcv] = generateBC_MFS_2D(ploton,boundunbound,corner1,corner2,nbc,b,Ubound,fos)

if nargin==7
    fos = 1;
end

% calculate center of mass of rectangle
com=[corner1(1)+corner2(1),corner1(2)+corner2(2)]/2;

switch boundunbound
    case 'unbounded'
        bcxcoor = linspace(corner1(1),corner2(1),nbc);
        bcycoor = linspace(corner1(2),corner1(2),nbc);
        % define speed of moving boundary
        bcu= linspace(0,0,nbc);
        bcv= linspace(0,0,nbc);
        velboundleft=max(round(nbc*(1/2-fos/2)),2);
        velboundright=min(round(nbc*(1/2+fos/2)),nbc);
        bcu(velboundleft-1:velboundright)=Ubound;
        % place stokeslets at a fractional distance b 
        stokescoordx=bcxcoor+b*(bcxcoor-com(1));
        stokescoordy=bcycoor+b*(bcycoor-com(2));
    case 'bounded'
        bxbot=linspace(corner1(1),corner2(1),nbc);
        bybot=linspace(corner1(2),corner1(2),nbc);
        bxright=linspace(corner2(1),corner2(1),nbc);
        byright=linspace(corner1(2),corner2(2),nbc);
        bxtop=linspace(corner2(1),corner1(1),nbc);
        bytop=linspace(corner2(2),corner2(2),nbc);
        bxleft=linspace(corner1(1),corner1(1),nbc);
        byleft=linspace(corner2(2),corner1(2),nbc);
        bxbot=bxbot(1:end-1);bybot=bybot(1:end-1);
        bxright=bxright(1:end-1);byright=byright(1:end-1);
        bxtop=bxtop(1:end-1);bytop=bytop(1:end-1);
        bxleft=bxleft(1:end-1);byleft=byleft(1:end-1);
        bcxcoor=[bxbot bxright bxtop bxleft];
        bcycoor=[bybot byright bytop byleft];
        
        % place stokeslets at a fractional distance b 
        stokescoordx=bcxcoor+b*(bcxcoor-com(1));
        stokescoordy=bcycoor+b*(bcycoor-com(2));
        
        % define speed of moving boundary
        bcv=zeros(size(stokescoordy));
        velboundleft=max(round(nbc*(1/2-fos/2)),2);
        velboundright=min(round(nbc*(1/2+fos/2)),nbc);
        bubot=linspace(0,0,nbc);
        bubot(velboundleft:velboundright)=Ubound;
        buright=linspace(0,0,nbc);butop=linspace(0,0,nbc);buleft=linspace(0,0,nbc);
        bubot=bubot(1:end-1);buright=buright(1:end-1);butop=butop(1:end-1);buleft=buleft(1:end-1);
        bcu=[bubot buright butop buleft];
end

% plot geometry if selected
if ploton==1
    figure(1);plot(stokescoordx,stokescoordy,'x',bcxcoor,bcycoor,'o',com(1),com(2),'+');
    figure(2);quiver(bcxcoor,bcycoor,bcu,bcv);
end

end