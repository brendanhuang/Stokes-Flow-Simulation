%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: calcUV_MFS_2D
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
% Description: This function estimates the flow field u and v
% by computing the forward equation u = G f in the method of fundamental
% solutions
% 
% Inputs: 
% qcoeff - strength / coefficients of stokeslets
% Amat - flow Green's function / susceptibility matrix
% xpts, ypts - number of xpts and y pts in the interior grid to calculate 
% Note this must be the same as the xpts, ypts used to calculate the
% original Amat
%
%
% Outputs:
% u,v - vectorial  flow field

function [u,v] = calcUV_MFS_2D(qcoeff,Amat,xpts,ypts)
uvtot=Amat*qcoeff; % Calculate forward equation
[xgrid,ygrid]=meshgrid(linspace(1,xpts,xpts),linspace(1,ypts,ypts));
u1=uvtot(1:xpts*ypts);
v1=uvtot(xpts*ypts+1:end);

u=accumarray([ygrid(:),xgrid(:)],u1,[ypts xpts]);
v=accumarray([ygrid(:),xgrid(:)],v1,[ypts xpts]);
end