%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: solve_mixedBC
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
% Description: This function solves the equation u = G f + T u
% where u is the flow on the boundary, f is the traction on the boundary, G
% is the single layer potential, and T is the double layer potential. Some
% elements of both u and f may be known, while others are unknown. 
% The solution is done by collecting all the known coeeficients in one vector and
% all the unknown coefficients in another vector (swapping traction and boundary conditions
% using the matlab deal function, inverting the equation, and then swaping the 
% solved for traction and boundary conditions back in. 
%
% Inputs: 
% BCtotal - known velocity / traction values on boundary, all dimensions (2D 
% or 3D) in a single column 
% alpha - index specifying which boundary element is traction versus velocity
% SLmat - single layer potential
% DLmat - double layer potential
% visc - viscosity 
% dimen - 2D vs. 3D 
%
% Outputs:
% u0 velocity values of boundary, all dimensions in a single column
% f0 traction values on boundaries, all dimensions in a single column
% bcu, bcv, [bcw] - u0 split into proper dimensions
% fx, fy, [fz] - f0 split into proper dimensions
% We have a variable out structure depending on the number of dimensions
% For 2D for example, bcu, bcv, fx, and fy are nout 1 - 4

function [u0,f0,nout1,nout2,nout3,nout4,nout5,nout6] = solve_mixedBC(SLmat,DLmat,BCtotal,visc,alpha,dimen)

% We start by multiplying our single and double layer potential matrices by
% the proper coefficients
switch dimen
    case '2D'
    A=1/(2*pi*visc)*SLmat;
    B=(1/(2*pi)*DLmat-eye(size(DLmat)));
    case '3D'
    A=1/(4*pi*visc)*SLmat;
    B=(1/(4*pi)*DLmat-eye(size(DLmat)));
end
% Our Boundary Condition vector BCtotal is a long vector that contains
% either values of u or f. Whether the value is u or f is determined by the
% alpha vector. Here, if a boundary condition is a traction one (f), then
% we swap the single and double layer potentials for that element
[A(:,alpha),B(:,alpha)]=deal(-B(:,alpha),-A(:,alpha));

% may want to A\B this in the future
f0=inv(A)*B*BCtotal; 
u0=BCtotal;

% f0 is our now solved for coefficients. Most are traction, but the
% elements in BC total that were traction need to be swapped back
[u0(alpha),f0(alpha)]=deal(f0(alpha),u0(alpha));

% Split u0 and f0 into respective x,y, [z] dimensions
switch dimen
    case '2D'
        n1=numel(u0)/2;
        nout1=u0(1:n1);
        nout2=u0(n1+1:end);
        nout3=f0(1:n1);
        nout4=f0(n1+1:end);
    case '3D'
        n1=numel(u0)/3;
        nout1=u0(1:n1);
        nout2=u0(n1+1:2*n1);
        nout3=u0(2*n1+1:end);
        nout4=f0(1:n1);
        nout5=f0(n1+1:2*n1);
        nout6=f0(2*n1+1:end);
end
end