%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: calc_Amat_MFS_2D
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
% Description: This function uses Green's functions in 2D to calculate
% the coupling between boundary elements and points. It does so for flow,
% pressure, and the streamfunction
% 
%
% Inputs:
% boundarycoordx,boundarycoordy - location of boundary points
% stokescoordx,stokescoordy - location of stokeslets
% visc - viscosity
% 
% Outputs: 
% Amat - Flow Green's function / susceptibility matrix
% Optional Outputs:
% Apres, Apsi - Pressure and streamfunction Green's functions /
% susceptibility matrix

function [Amat,Apres,Apsi] = calc_Amat_MFS_2D(boundarycoordx,boundarycoordy,stokescoordx,stokescoordy,visc)
if nargin==4
    visc=1;
end


nstokes=numel(stokescoordx);
npts=numel(boundarycoordx);

% calculate displacement vector and distances between stokselets / points
xij=transpose(repmat(boundarycoordx,nstokes,1))-repmat(stokescoordx,npts,1);
yij=transpose(repmat(boundarycoordy,nstokes,1))-repmat(stokescoordy,npts,1);
rij=(xij.^2+yij.^2).^0.5;

% calculate the xx, xy, yx, yy components of the matrix
fij=-2*log(rij)+2*(xij.^2).*(rij.^(-2))-3;
gij=2*xij.*yij.*(rij.^(-2));
hij=2*xij.*yij.*(rij.^(-2));
kij=-2*log(rij)+2*(yij.^2).*(rij.^(-2))-3;

% put all components back together. Note here, our Amat is slightly
% different from our BEM matrix because it contains the factor 8 pi visc
% already incorporated.  Thus when we solve the inverse equation, we don't
% need to weight by those factors.
Amat=[fij, gij; hij, kij]/(8*pi*visc);

% calculate the pressure Green's function matrix
if nargout>1
    fpres=xij.*rij.^(-2);
    gpres=yij.*rij.^(-2);
    Apres=[fpres, gpres]/(2*pi);
end

% calculate the streamfunction Green's function matrix
if nargout>2
   fpsi=-2*yij.*log(rij)-yij;
   gpsi=2*xij.*log(rij)-xij;
   Apsi=[fpsi,gpsi]/(8*pi*visc);
end

end
