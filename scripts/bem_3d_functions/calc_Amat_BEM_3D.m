%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: calc_Amat_BEM_3D
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
% Description: This function uses Green's functions in 3D to calculate
% the coupling between boundary elements and points. This function simply
% loops the calc_SLDL_3D function, a function that calculates individual
% elements in the matrix.
% 
%
% Inputs:
% vertpts - location of vertices of triangles of boundary elements 
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
% Both have dimensions of 3m × 3n

function [SLmat,DLmat] =calc_Amat_BEM_3D(p,vertpts,elemvert,normv,p_on)
% Pre-allocate memory for a 9 different m x n matrices, and then we will
% put them back together at the final step

n=size(elemvert,1);
m=size(p,1);

SLxx=zeros(m,n);
SLxy=zeros(m,n);
SLxz=zeros(m,n);
SLyx=zeros(m,n);
SLyy=zeros(m,n);
SLyz=zeros(m,n);
SLzx=zeros(m,n);
SLzy=zeros(m,n);
SLzz=zeros(m,n);

DLxx=zeros(m,n);
DLxy=zeros(m,n);
DLxz=zeros(m,n);
DLyx=zeros(m,n);
DLyy=zeros(m,n);
DLyz=zeros(m,n);
DLzx=zeros(m,n);
DLzy=zeros(m,n);
DLzz=zeros(m,n);
    for i=1:m 
        for j=1:n           
            % To evaluate the Green's function between each boundary
            % element and point, we need the vertices of the boundary
            % element (r1, r2, r3), the normal vector of the boundary element
            % (normal), and the coordinate of the point (p1)
            r1=vertpts(elemvert(j,1),:);
            r2=vertpts(elemvert(j,2),:);
            r3=vertpts(elemvert(j,3),:);
            normal=normv(j,:);
            p1=p(i,:);            
            [GSL,GDL]=calc_SLDL_3D(p1,r1,r2,r3,normal,p_on&i==j);
            
            SLxx(i,j)=GSL(1,1);
            SLxy(i,j)=GSL(1,2);
            SLxz(i,j)=GSL(1,3);
            SLyx(i,j)=GSL(2,1);
            SLyy(i,j)=GSL(2,2);
            SLyz(i,j)=GSL(2,3);
            SLzx(i,j)=GSL(3,1);
            SLzy(i,j)=GSL(3,2);
            SLzz(i,j)=GSL(3,3);
            
            DLxx(i,j)=GDL(1,1);
            DLxy(i,j)=GDL(1,2);
            DLxz(i,j)=GDL(1,3);
            DLyx(i,j)=GDL(2,1);
            DLyy(i,j)=GDL(2,2);
            DLyz(i,j)=GDL(2,3);
            DLzx(i,j)=GDL(3,1);
            DLzy(i,j)=GDL(3,2);
            DLzz(i,j)=GDL(3,3);
        end
    end
    
SLmat=[SLxx,SLxy,SLxz;SLyx,SLyy,SLyz;SLzx,SLzy,SLzz];
DLmat=[DLxx,DLxy,DLxz;DLyx,DLyy,DLyz;DLzx,DLzy,DLzz];

end