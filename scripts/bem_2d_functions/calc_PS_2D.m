%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: calc_PS_2D
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
% Description: This function uses the pressure and shear stress
% Green's functions in 2D to calculate the coupling between boundary 
% elements and points. It fills in the elements of the full Single Layer
% and Double Layer potential matrices.
% 
%
% Inputs:
% p1 point to calculate 
% r1, r2 vertices of line segment (2D) of boundary element; 
% normv normal vectors of boundary element
% p_on is our point being evaluated on our boundary element?
% 
% Outputs: 
% PSL, PDL: single and double layer potential coupling for pressure
% SSL, SDL: single and double layer potential coupling for shear stress. 
%
% References: for form of Green's Functions, see
% [1] C. Pozrikidis. A practical guide to boundary element methods 
% with the software library BEMLIB. Chapman & Hall/CRC, Boca Raton, 2002.
% Ch. 7


function [PSL,PDL,SSL,SDL] = calc_PS_2D(p1,r1,r2,normal)

% To evaluate the Green's function on a single element, we use a
% Gauss-Legendre Quadrature. gaussleg, defined below, gives the weights and evaluation locations of
% the Gauss-Legendre quadrature method
[w,x,~]=gaussleg;

% w and x come out assuming a line segment from 0 - 1. We need to reshift coordinates 
xx=r1(1)+x*(r2(1)-r1(1));
yy=r1(2)+x*(r2(2)-r1(2));
A=norm(r2-r1);

xw12=repmat(w,[2 1]);
xw222=permute(repmat(w',1,2,2,2),[2,3,4,1]);

% pressure_single_2d and pressure_double_2d  evaluates the single 
% and double layer potential at the series of points
% xx, yy on the boundary element, and coupled to the points p1.
% lponq is no longer necessary because we don't have the analytical
% expressions to deal with singularities
PSL = pressure_single_2d(xx,yy,p1(1),p1(2));
PSL = A*dot(PSL,xw12,2);
PDL = pressure_double_2d(xx,yy,p1(1),p1(2),normal);
PDL = A*dot(PDL,xw12,2);

% only calculate shear if we need to
if nargout>2
SSL = stress_single_2d(xx,yy,p1(1),p1(2));
SSL = A*dot(SSL,xw222,4);
SDL = stress_double_2d(xx,yy,p1(1),p1(2),normal);
SDL = A*dot(SDL,xw222,4);
end

 
end

function [wts,pts,n]=gaussleg() 
% gaussleg gives Gauss-Legendre Quadrature 8 points and weights to approximate 
% an integral
n=8;
wts= [ 5.061426814519E-02 0.111190517227 0.156853322939 0.181341891689 0.181341891689 0.156853322939 0.111190517227 5.061426814519E-02];
pts= [ 1.985507175123E-02 0.101666761293 0.237233795042 0.408282678752 0.591717321248 0.762766204958 0.898333238707 0.980144928249];
end

function PS2D=pressure_single_2d(x,y,x0,y0)
% Here we use the pressure Green's function for points in our Gaussian
% Quadrature we calculate the displacement vector dx, dy, and the 
% distance dr, and feed it into the formula
% See Ref [1], eq 7.2.13
      dx = x-x0;
      dy = y-y0;
      dr = sqrt(dx.^2+dy.^2);
      PS2D(1,:) = 2*dx./dr.^2;
      PS2D(2,:) = 2*dy./dr.^2;
end

function PD2D=pressure_double_2d(x,y,x0,y0,normal)
% Here we do the same for the double layer potential for pressure
% See Ref [1], eq 7.3.8
      dx = x-x0;
      dy = y-y0;
      dr = sqrt(dx.^2+dy.^2);
      np = numel(dx);
      PD=zeros(2,2,np);
      PD2D=zeros(2,np);
      PD(1,1,:)=2*(-1./dr.^2+2*dx.^2./dr.^4);
      PD(1,2,:)=2*(2*dx.*dy./dr.^4);
      PD(2,1,:)=2*(2*dx.*dy./dr.^4);
      PD(2,2,:)=2*(-1./dr.^2+2*dy.^2./dr.^4);

      % we have to dot with our normal vector again
      for j=1:2
            for k=1:np
            PD2D(j,k)=normal*squeeze(PD(j,:,k))';
            end
      end
end

 function Tijk=stress_single_2d(x,y,x0,y0)
% Here we plug in the Green's function for our shear stress to
% get our single layer potential. Because of the numerous components, we
% chose again to loop through the dimensions
% See Ref [1], eq. 7.2.13
      dx = x-x0;
      dy = y-y0;
      dr = sqrt(dx.^2+dy.^2);
      np = numel(dx);
      dt=[dx;dy];
     Tijk=zeros(2,2,2,np);
     for i=1:2
         for j=1:2
             for k=1:2
                 Tijk(i,j,k,:)=-4*dt(i,:).*dt(j,:).*dt(k,:)./dr.^4;
             end
         end
     end
end

function Kijk=stress_double_2d(x,y,x0,y0,normal)
% Here we plug in the Green's function for our shear stress to
% get our double layer potential. Because of the numerous components, we
% chose again to loop through the dimensions
% See Ref [1], Eq. 7.3.11
      dx = x-x0;
      dy = y-y0;
      np=numel(dx);         
      dr=sqrt(dx.^2+dy.^2);
      dt=[dx;dy];

    K_ijkl=zeros(2,2,2,2,np);
    Kijk=zeros(2,2,2,np);

    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:2
K_ijkl(i,j,k,l,:)=KD(i,k).*KD(j,l).*4./dr.^2+4./dr.^4.*(KD(i,j).*dt(k,:).*dt(l,:)+ ...
    KD(i,l).*dt(k,:).*dt(j,:)+KD(k,j).*dt(i,:).*dt(l,:)+KD(k,l).*dt(i,:).*dt(j,:)) ...
    -32.*dt(i,:).*dt(j,:).*dt(k,:).*dt(l,:)./dr.^6;
                end
            end
        end
    end

% we will contract over our normal vector again    
    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:np
            Kijk(i,j,k,l)=normal*squeeze(K_ijkl(i,j,k,:,l));
                end
            end
        end
    end
end

function k = KD(a,b)
% Here we are defining our Kronecker delta function (KD). It is employed in
% calculating the shear stress potential
    if a == b
        k = 1;
    else
        k = 0;
    end
end
