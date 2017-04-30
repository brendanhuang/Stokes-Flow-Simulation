%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: calc_Amat_BEM_2D
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
% the coupling between boundary elements and points. This function simply
% loops the calc_SLDL_2D function, a function that calculates individual
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
% SLmat single layer potential susceptibility matrix; 
% DL- mat double layer potential susceptibility matrix.
% Both have dimensions of 2m × 2n


function [SLmat,DLmat] = calc_Amat_BEM_2D(p,vertpts,elemvert,normv,p_on)
% Pre-allocate memory for a 2m x 2n matrix
m=size(p,1);
n=size(elemvert,1);
SLxx=zeros(m,n);
SLxy=zeros(m,n);
SLyx=zeros(m,n);
SLyy=zeros(m,n);
DLxx=zeros(m,n);
DLxy=zeros(m,n);
DLyx=zeros(m,n);
DLyy=zeros(m,n);

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
            [SLxx(i,j),SLxy(i,j),SLyx(i,j),SLyy(i,j),DLxx(i,j),DLxy(i,j),DLyx(i,j),DLyy(i,j)]=calc_SLDL_2D(p1,r1,r2,normal,p_on&i==j);
        end
    end
SLmat=[SLxx,SLxy;SLyx,SLyy];
DLmat=[DLxx,DLxy;DLyx,DLyy];
end