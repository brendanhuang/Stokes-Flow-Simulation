%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: calcPPsi_MFS_2D
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
% Description: This function estimates either the streamfunction or
% pressure field by computing the forward equation in MFS
% 
% Inputs: 
% qcoeff - strength / coefficients of stokeslets
% Pint - flow Green's function / susceptibility matrix for either pressure
% or streamfunction
% xpts, ypts - number of xpts and y pts in the interior grid to calculate 
% Note this must be the same as the xpts, ypts used to calculate the
% original Amat
%
%
% Outputs:
% p - field for either pressure or streamfunction

function p = calcPPsi_MFS_2D(qcoeff,Pint,xpts,ypts)

ptot=Pint*qcoeff;
[xgrid,ygrid]=meshgrid(linspace(1,xpts,xpts),linspace(1,ypts,ypts));

p=accumarray([ygrid(:),xgrid(:)],ptot,[ypts xpts]);

end