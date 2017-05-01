%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: grid2tri
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
% Description: This function takes rectangular grid defined by the number
% of points along each dimension (nx = s1, ny = s2) rectangle into two right triangles 
% (upper right (ur) and bottom left (bl)). It attempts to orient both ur and
% bl triangles such that the surfaces both have normals oriented in the same 
% direction.
% 
% Inputs: 
% s1, s2 - number of x and y points in each dimension
% optional: flipnorm - in case the normal vector is oriented in the wrong direction
% (i.e. facing out of the cube), flipnorm reverses the direction of the
% unit normal vector
%
% Outputs:
% tribl, triur: bottom left and upper right triangles consisting of the
% indices to take from a mesh grid.

function [tribl,triur] = grid2tri(s1,s2,flipnorm)
if nargin==2
    flipnorm=1;
end


tribl=zeros((s1-1)*(s2-1),3);
triur=zeros((s1-1)*(s2-1),3);
i3=1;
for i1=1:s1-1
    for i2=1:s2-1
    tribl(i3,1)=sub2ind([s1,s2],i1,i2);
    tribl(i3,2)=sub2ind([s1,s2],i1+1,i2);
    tribl(i3,3)=sub2ind([s1,s2],i1,i2+1);
    triur(i3,2)=sub2ind([s1,s2],i1+1,i2);
    triur(i3,1)=sub2ind([s1,s2],i1,i2+1);
    triur(i3,3)=sub2ind([s1,s2],i1+1,i2+1);
    i3=i3+1;
    end
end
if flipnorm==2
    triur=[triur(:,2),triur(:,1),triur(:,3)];
    tribl=[tribl(:,2),tribl(:,1),tribl(:,3)];
end

end
