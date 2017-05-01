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

function [PSLmat,PDLmat,SSLmat,SDLmat] =calc_PSmat_BEM_3D(p,vertpts,elemvert,normv,p_on)

n=size(elemvert,1);
m=size(p,1);

% if we want to evaluate pressure and stress on the boundary, we need to
% deal with singularities. These cases have not yet been solved, so we
% return a matrix of zeros
if p_on
    PSLmat=zeros(m,3*n);
    PDLmat=zeros(m,3*n);
    SSLmat=zeros(9*m,3*n);
    SDLmat=zeros(9*m,3*n);
    disp('Not valid at boundaries');
elseif nargout<3
    % If we only want pressure and not shear stress, we'll only calculate
    % pressure

    PSLmat=zeros(3,m,n);
    PDLmat=zeros(3,m,n);
    SSLmat=zeros(9*m,3*n);
    SDLmat=zeros(9*m,3*n);

    for i=1:m 
        for j=1:n           
            % To evaluate the Green's function between each boundary
            % element and point, we need the endpoints of the boundary
            % element (r1, r2), the normal vector of the boundary element
            % (normal), and the coordinate of the point (p1)

            r1=vertpts(elemvert(j,1),:);
            r2=vertpts(elemvert(j,2),:);
            r3=vertpts(elemvert(j,3),:);
            normal=normv(j,:);
            p1=p(i,:);
            [PSLmat(:,i,j),PDLmat(:,i,j)] = calc_PS_3D(p1,r1,r2,r3,normal);
        end
    end
PSLmat=reshape(permute(PSLmat,[2,3,1]),[m,3*n]);
PDLmat=reshape(permute(PDLmat,[2,3,1]),[m,3*n]);
else
    
    % The default case is to calculate pressure and stress. In this case, to
    % make sure we have our indexing right, we will create a 3 and 5
    % dimensional matrix, but then reshape at the end to m x 3n and 9m x 3n 

    PSLmat=zeros(3,m,n);
    PDLmat=zeros(3,m,n);
    SSLmat=zeros(3,3,3,m,n);
    SDLmat=zeros(3,3,3,m,n);

    for i=1:m 
        for j=1:n           
            % To evaluate the Green's function between each boundary
            % element and point, we need the endpoints of the boundary
            % element (r1, r2), the normal vector of the boundary element
            % (normal), and the coordinate of the point (p1)
            r1=vertpts(elemvert(j,1),:);
            r2=vertpts(elemvert(j,2),:);
            r3=vertpts(elemvert(j,3),:);
            normal=normv(j,:);
            p1=p(i,:);            
            [PSLmat(:,i,j),PDLmat(:,i,j),SSLmat(:,:,:,i,j),SDLmat(:,:,:,i,j)] = calc_PS_3D(p1,r1,r2,r3,normal);            
        end
    end
PSLmat=reshape(permute(PSLmat,[2,3,1]),[m,3*n]);
PDLmat=reshape(permute(PDLmat,[2,3,1]),[m,3*n]);
SSLmat=reshape(permute(SSLmat,[4,1,2,5,3]),[9*m,3*n]);
SDLmat=reshape(permute(SDLmat,[4,1,2,5,3]),[9*m,3*n]);
 
end

end