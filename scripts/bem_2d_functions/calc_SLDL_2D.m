%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: calc_SLDL_2D
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
% the coupling between boundary elements and points. It fills in the
% elements of the full Single Layer and Double Layer potential matrices
% 
%
% Inputs:
% p1 point to calculate 
% r1, r2 vertices of line segment (2D) of boundary element; 
% normv normal vectors of boundary element
% p_on is our point being evaluated on our boundary element?
% 
% Outputs: 
% SLxx, SLxy, SLyx, SLyy: single layer potential coupling coefficients
% DLxx, DLxy, DLyx, DLyy double layer potential coupling coeffcients. 
% For both, there are 4 components
%
% References: for form of Green's Functions, see
% [1] C. Pozrikidis. A practical guide to boundary element methods 
% with the software library BEMLIB. Chapman & Hall/CRC, Boca Raton, 2002,
% Chapter 7


function [SLxx,SLxy,SLyx,SLyy,DLxx,DLxy,DLyx,DLyy] = calc_SLDL_2D(p1,r1,r2,normal,lponq)
% To evaluate the Green's function on a single element, we use a
% Gauss-Legendre Quadrature. gaussleg, defined below, gives the weights and evaluation locations of
% the Gauss-Legendre quadrature method
[w,x,~]=gaussleg;

% w and x come out assuming a line segment from 0 - 1. We need to reshift coordinates 
xx=r1(1)+x*(r2(1)-r1(1));
yy=r1(2)+x*(r2(2)-r1(2));
A=norm(r2-r1);

% gf_2d_fs now evaluates the single layer potential at the series of points
% xx, yy on the boundary element, and coupled to the 
% lponq deals with whether the point is collocalized on the boundary
% element itself, in which case we have to deal with a singularity
[Gxx,Gxy,Gyx,Gyy]=gf_2d_fs(xx,yy,p1(1),p1(2),lponq);

SLxx=A*(w*Gxx');
SLxy=A*(w*Gxy');
SLyx=A*(w*Gyx');
SLyy=A*(w*Gyy');

% In this 2D implementation, evaluating the Singularity done as per Pozrikidis
% When we evaluate gf_2d_fs, we simply remove the log(r) term. We then 
% add in an analytic evaluation now
if lponq
    a=norm(p1-r1);
    b=norm(p1-r2);
    SLxx=SLxx+(a+b-(a*log(a)+b*log(b))); 
    SLyy=SLyy+(a+b-(a*log(a)+b*log(b)));
end

% Now we calculate double layer potential.  Fortunately, the double layer potential
% is zero if we have a singularity because the dot product of the normal vector with
% the displacement vector between the point and the boundary element is zero (they are
% orthogonal or (n,u) = 0)
% If they are not collocalized, then we evaluate the stresslet Green's function
% at our specified Gauss-Legendre points
if lponq
DLxx=0;
DLxy=0;
DLyx=0;
DLyy=0;
else
[Sxx,Sxy,Syx,Syy]=stresslet_2d(xx,yy,p1(1),p1(2),normal);
 DLxx=A*(w*Sxx);
 DLxy=A*(w*Sxy);
 DLyx=A*(w*Syx);
 DLyy=A*(w*Syy);
end    

end

function [wts,pts,n]=gaussleg() 
% gaussleg gives Gauss-Legendre Quadrature 8 points and weights to approximate 
% an integral
n=8;
wts= [ 5.061426814519E-02 0.111190517227 0.156853322939 0.181341891689 0.181341891689 0.156853322939 0.111190517227 5.061426814519E-02];
pts= [ 1.985507175123E-02 0.101666761293 0.237233795042 0.408282678752 0.591717321248 0.762766204958 0.898333238707 0.980144928249];
end

function [Gxx,Gxy,Gyx,Gyy]=gf_2d_fs(x,y,x0,y0,lponq)
% Here we use the Oseen-Burgers tensor for points in our Gaussian
% Quadrature
% See Ref [1], Eq. 7.2.13

% Calculate the displacement vector dx, dy, and the distance dr
      dx = x-x0;
      dy = y-y0;
      dr = sqrt(dx.^2+dy.^2);

% As noted above, if the point is in the boundary element, we ignore the
% log(r) term for now. We'll add back in the analytic answer later
      if lponq
      Gxx =dx.^2./dr.^2;
      Gxy =dx.*dy./dr.^2;
      Gyx =dx.*dy./dr.^2;
      Gyy =dy.^2./dr.^2;
      else
      Gxx = -log(dr) + dx.^2./dr.^2;
      Gxy = dx.*dy./dr.^2;
      Gyx = dx.*dy./dr.^2;
      Gyy = -log(dr) + dy.^2./dr.^2;
      end
end

function [Txx,Txy,Tyx,Tyy]=stresslet_2d(x,y,x0,y0,normal)
% Here we plug in the Green's function for our stresslet to
% get our double layer potential. 
% See Ref [1], Eq. 7.2.13
      dx = x-x0;
      dy = y-y0;
      np=numel(dx);
         
      dr=sqrt(dx.^2+dy.^2);
      dt=[dx;dy];


    Tijk=zeros(2,2,2,np);
    Tij=zeros(2,2,np);

    % There are a number of components to evaluate, so in order to keep the
    % scripting as compact as possible and to avoid scripting errors, we
    % use a general expression and loop through the indices
    for i=1:2
        for j=1:2
            for k=1:2
                Tijk(i,j,k,:)=-4*dt(i,:).*dt(j,:).*dt(k,:)./dr.^4;
            end
        end
    end

    % Now we take the inner product of Tijk with nk, our normal vector
    for i=1:2
        for j=1:2
            for k=1:np
            Tij(i,j,k)=normal*squeeze(Tijk(i,j,:,k));
            end
        end
    end
    Txx=squeeze(Tij(1,1,:));
    Txy=squeeze(Tij(1,2,:));
    Tyx=squeeze(Tij(2,1,:));
    Tyy=squeeze(Tij(2,2,:));
end