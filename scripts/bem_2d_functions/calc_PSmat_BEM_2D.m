%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: calc_PSmat_BEM_2D
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
% Description: This function uses Green's functions for pressure and 
% shear stress in 2D to calculate the coupling between boundary elements and points. This function simply
% loops the calc_PSmat_BEM_2D function, a function that calculates individual
% elements in the matrix.
% 
%
% Inputs:
% vertpts - location of vertices of line (2D) of boundary elements 
% elemvert - vector to identify which vertices belong to which boundary element
% p - location of points that are coupled to boundary elements are evaluated
% p_on - a true / false variable letting us know if we expect some of the
% points that we evluate to be collocated with boundary elements, thus
% giving rise to a possible singularity
% normv - normal vector of each boundary element
% 
% Outputs: 
% PSLmat - single layer potential pressure matrix 
% PDLmat - double layer potential pressure matrix
% SSLmat - single layer potential shear stress matrix 
% SDLmat - double layer potential shear stress matrix
% The pressure values are scalar, while the shear stress is a d x d tensor
% Thus, the matrices are m x 2n and 4 m x 2n for d=2.

function [PSLmat,PDLmat,SSLmat,SDLmat] = calc_PSmat_BEM_2D(p,vertpts,elemvert,normv,p_on)
m=size(p,1);
n=size(elemvert,1);

% if we want to evaluate pressure and stress on the boundary, we need to
% deal with singularities. These cases have not yet been solved, so we
% return a matrix of zeros
if p_on
    PSLmat=zeros(m,2*n);
    PDLmat=zeros(m,2*n);
    SSLmat=zeros(4*m,2*n);
    SDLmat=zeros(4*m,2*n);
    
% If we only want pressure and not shear stress, we'll only calculate
% pressure
elseif nargout<3
    PSLmat=zeros(2,m,n);
    PDLmat=zeros(2,m,n);
    SSLmat=zeros(4*m,2*n);
    SDLmat=zeros(4*m,2*n);
    
    for i=1:m 
        for j=1:n
            r1=vertpts(elemvert(j,1),:);
            r2=vertpts(elemvert(j,2),:);
            normal=normv(j,:);
            p1=p(i,:);            
            [PSLmat(:,i,j),PDLmat(:,i,j),SSLmat(:,:,:,i,j),SDLmat(:,:,:,i,j)]=calc_PS_2D(p1,r1,r2,normal,p_on&i==j);
        end
    end
    PSLmat=reshape(permute(PSLmat,[2,3,1]),[m,2*n]);
    PDLmat=reshape(permute(PDLmat,[2,3,1]),[m,2*n]);

else
% The default case is to calculate pressure and stress. In this case, to
% make sure we have our indexing right, we will create a 3 and 5
% dimensional matrix, but then reshape at the end to m x 2n and 4m x 2n 
    PSLmat=zeros(2,m,n);
    PDLmat=zeros(2,m,n);
    SSLmat=zeros(2,2,2,m,n);
    SDLmat=zeros(2,2,2,m,n);
    for i=1:m 
        for j=1:n
            % To evaluate the Green's function between each boundary
            % element and point, we need the endpoints of the boundary
            % element (r1, r2), the normal vector of the boundary element
            % (normal), and the coordinate of the point (p1)

            r1=vertpts(elemvert(j,1),:);
            r2=vertpts(elemvert(j,2),:);
            normal=normv(j,:);
            p1=p(i,:);            
            [PSLmat(:,i,j),PDLmat(:,i,j),SSLmat(:,:,:,i,j),SDLmat(:,:,:,i,j)]=calc_PS_2D(p1,r1,r2,normal);
        end
    end
PSLmat=reshape(permute(PSLmat,[2,3,1]),[m,2*n]);
PDLmat=reshape(permute(PDLmat,[2,3,1]),[m,2*n]);
SSLmat=reshape(permute(SSLmat,[4,1,2,5,3]),[4*m,2*n]);
SDLmat=reshape(permute(SDLmat,[4,1,2,5,3]),[4*m,2*n]);
end

end