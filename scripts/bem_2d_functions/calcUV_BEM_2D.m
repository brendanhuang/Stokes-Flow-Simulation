%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: calcUV_BEM_2D
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
% by computing the forward equation u_int = G f_boundary + T u_boundary,
% where u_int is the flow field inside the domain, G is the single layer
% potential, T is the double layer potential, f_boundary is the traction at
% the boundary, and u_boundary is the flow at the boundary. It can also
% calculate the analgous pressure field (press) and shear stress tensor
% (stress) using the single and doulbe layer potentials
% 
% Inputs: 
% f0, u0 - traction and velocity on boundaries
% SLmat, DLmat - single and double layer potential susceptibility matrices
% linking elements on boundary with interior points.
% visc - viscosity
% xpts, ypts: number of points to calculate in the interior. Note this must
% be the same as the xpts, ypts used to calculate the single and double
% layer potentials in calc_Amat_BEM_2D. The position of the points is
% implicit in those matrices.
% Optional Inputs: 
% PSL, PDL, SSL, SDL single and double layer potentials for pressure and stress, 
% if we would like to directly calculate these quantities on the interior
%
% Outputs:
% u,v - vectorial  flow field
% Optional Outputs: 
% pres, stress - presssure field, shear stress tensor

function [u,v,pres,stress] = calcUV_BEM_2D(f0,u0,SLmat,DLmat,xpts,ypts,visc,PSL,PDL,SSL,SDL)
    
% calculate the forward equation 
uvtot = -1/(4*pi*visc)*SLmat*f0+1/(4*pi)*DLmat*u0;

% use the information of our number of xpts and ypts to reshape our flow
% field properly
[xgrid,ygrid]=ndgrid(linspace(1,xpts,xpts),linspace(1,ypts,ypts));
u1=uvtot(1:xpts*ypts);
v1=uvtot(xpts*ypts+1:end);

u=accumarray([xgrid(:),ygrid(:)],u1,[xpts ypts]);
v=accumarray([xgrid(:),ygrid(:)],v1,[xpts ypts]);

% if desired, calculate the forward equation for pressure
if nargout==3
    pres1 = -1/(4*pi)*PSL*f0+visc/(4*pi)*PDL*u0;
    pres=accumarray([xgrid(:),ygrid(:)],pres1,[xpts ypts]); 

% if desired, calculate the forward equation for shear stress as well
elseif nargout ==4
    pres1 = -1/(4*pi)*PSL*f0+visc/(4*pi)*PDL*u0;
    stress1 = -1/(4*pi)*SSL*f0+visc/(4*pi)*SDL*u0;
    pres=accumarray([xgrid(:),ygrid(:)],pres1,[xpts ypts]); 
    stress1=reshape(stress1,[xpts*ypts,2,2]);
    stress=zeros(xpts,ypts,2,2);

%reshape stress into tensor of xpts x ypts x 2 x 2
    for i=1:2
       for j=1:2
           stress(:,:,i,j)=accumarray([xgrid(:),ygrid(:)],squeeze(stress1(:,i,j)),[xpts ypts]);
       end
    end
end

end