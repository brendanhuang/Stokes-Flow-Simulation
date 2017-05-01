%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: calc_PS_3D
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
% Green's functions in 3D to calculate the coupling between boundary 
% elements and points. It fills in the elements of the full Single Layer
% and Double Layer potential matrices.
% 
%
% Inputs:
% p1 point to calculate 
% r1, r2, r3 - vertices of triangle of boundary element; 
% normv - normal vectors of boundary element
% p_on - is our point being evaluated on our boundary element?
% 
% Outputs: 
% PSL, PDL: single and double layer potential coupling for pressure
% SSL, SDL: single and double layer potential coupling for shear stress. 
%
% References: for form of Green's Functions, see
% [1] C. Pozrikidis. A practical guide to boundary element methods 
% with the software library BEMLIB. Chapman & Hall/CRC, Boca Raton, 2002.
% Ch. 7

function [PSL,PDL,SSL,SDL] = calc_PS_3D(p1,r1,r2,r3,normal)
% We need to map a 2D triangle that lives in 3 dimensions onto a a flat 2D
% coordinate system before we do our Gaussian quadrature. We take advantage
% of the fact that we are in a cube by first finding out which dimension
% has the same coordinates for all points, then define "x" and "y" accordingly
whichone=whichissame(r1,r2,r3);
% Calculate the area A of the triangle
A=norm(cross(r2-r1,r3-r1))/2.0;

% To evaluate the Green's function on a single element, we use a
% Gauss-Legendre Quadrature. gaussleg, defined below, gives the weights and evaluation locations of
% the Gauss-Legendre quadrature method


    N=5;
    xw = TriGL(N);
    if whichone==1
        x = r1(1)*ones(size(xw(:,1)));
        y = r1(2)*(1-xw(:,1)-xw(:,2))+r2(2)*xw(:,1)+r3(2)*xw(:,2);
        z = r1(3)*(1-xw(:,1)-xw(:,2))+r2(3)*xw(:,1)+r3(3)*xw(:,2);
    elseif whichone==2
        x = r1(1)*(1-xw(:,1)-xw(:,2))+r2(1)*xw(:,1)+r3(1)*xw(:,2);
        y = r1(2)*ones(size(xw(:,1)));
        z = r1(3)*(1-xw(:,1)-xw(:,2))+r2(3)*xw(:,1)+r3(3)*xw(:,2);
    elseif whichone==3
        x = r1(1)*(1-xw(:,1)-xw(:,2))+r2(1)*xw(:,1)+r3(1)*xw(:,2);
        y = r1(2)*(1-xw(:,1)-xw(:,2))+r2(2)*xw(:,1)+r3(2)*xw(:,2);
        z = r1(3)*ones(size(xw(:,1)));
    end
    x=x';y=y';z=z';
    xw13=repmat(xw(:,3)',[3 1]);
    xw333=permute(repmat(xw(:,3),1,3,3,3),[2,3,4,1]);

% pressure_single_3d and pressure_double_3d  evaluate the single 
% and double layer potential at the series of points
% xx, yy on the boundary element, and coupled to the points p1.
% lponq is no longer necessary because we don't have the analytical
% expressions to deal with singularities

    PSL = pressure_single_3d(x,y,z,p1(1),p1(2),p1(3));
    PSL = A*dot(PSL,xw13,2);
    PDL = pressure_double_3d(x,y,z,p1(1),p1(2),p1(3),normal);
    PDL = A*dot(PDL,xw13,2);

% only calculate shear if we need to
    if nargout>2
    SSL = stress_single_3d(x,y,z,p1(1),p1(2),p1(3));
    SSL = A*dot(SSL,xw333,4);
    SDL = stress_double_3d(x,y,z,p1(1),p1(2),p1(3),normal);
    SDL = A*dot(SDL,xw333,4);
    end
    
end

function whichone=whichissame(r1,r2,r3)
% simple function to determine which of the three coordinates is invariant
% across all the triangles
if r1(1)==r2(1) && r1(1)==r3(1)
    whichone=1;
elseif r1(2)==r2(2) && r1(2)==r3(2)
    whichone=2;
elseif r1(3)==r2(3) && r1(3)==r3(3)
    whichone=3;
end
end

function PS3D=pressure_single_3d(x,y,z,x0,y0,z0)
% Here we use the pressure Green's function for points in our Gaussian
% Quadrature we calculate the displacement vector dx, dy, and the 
% distance dr, and feed it into the formula
% See Ref [1], eq 7.2.14
      dx = x-x0;
      dy = y-y0;
      dz = z-z0;
      dr = sqrt(dx.^2+dy.^2+dz.^2);
      np = numel(dx);
      PS3D=zeros(3,np);
      PS3D(1,:)=2*dx./dr.^3;
      PS3D(2,:)=2*dy./dr.^3;
      PS3D(3,:)=2*dz./dr.^3;
end

function PD3D=pressure_double_3d(x,y,z,x0,y0,z0,normal)
% Here we do the same for the double layer potential for pressure
% See Ref [1], eq 7.4.4
      dx = x-x0;
      dy = y-y0;
      dz = z-z0;
      dr = sqrt(dx.^2+dy.^2+dz.^2);
      np = numel(dx);
      PD=zeros(3,3,np);
      PD3D=zeros(3,np);
      PD(1,1,:)=4*(-1./dr.^3+3*dx.^2./dr.^5);
      PD(1,2,:)=4*(3*dx.*dy./dr.^5);
      PD(1,3,:)=4*(3*dx.*dz./dr.^5);
      PD(2,1,:)=4*(3*dx.*dy./dr.^5);
      PD(2,2,:)=4*(-1./dr.^3+3*dy.^2./dr.^5);
      PD(2,3,:)=4*(3*dy.*dz./dr.^5);
      PD(3,1,:)=4*(3*dx.*dz./dr.^5);
      PD(3,2,:)=4*(3*dx.*dy./dr.^5);
      PD(3,3,:)=4*(-1./dr.^3+3*dz.^2./dr.^5);
        for j=1:3
            for k=1:np
                %dot with normal vector
            PD3D(j,k)=normal*squeeze(PD(j,:,k))';
            end
        end
end

 function Tijk=stress_single_3d(x,y,z,x0,y0,z0)
 % Here we plug in the Green's function for our shear stress to
% get our single layer potential. Because of the numerous components, we
% chose again to loop through the dimensions
% See Ref [1], eq. 7.2.14

      dx = x-x0;
      dy = y-y0;
      dz = z-z0;
      dr = sqrt(dx.^2+dy.^2+dz.^2);
      np = numel(dx);
      dt=[dx;dy;dz];
     Tijk=zeros(3,3,3,np);
     for i=1:3
         for j=1:3
             for k=1:3
                 Tijk(i,j,k,:)=-6*dt(i,:).*dt(j,:).*dt(k,:)./dr.^5;
             end
         end
     end
end

function Kijk=stress_double_3d(x,y,z,x0,y0,z0,normal)
% Here we plug in the Green's function for our shear stress to
% get our double layer potential. Because of the numerous components, we
% chose again to loop through the dimensions
% See Ref [1], Eq. 7.4.6

      dx = x-x0;
      dy = y-y0;
      dz = z-z0;
      np=numel(dx);         
      dr=sqrt(dx.^2+dy.^2+dz.^2);
      dt=[dx;dy;dz];


    Psi_ijkl=zeros(3,3,3,3,np);
    Kijk=zeros(3,3,3,np);

    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
Psi_ijkl(i,j,k,l,:)=KD(i,k).*KD(j,l).*4./dr.^3+6./dr.^5.*(KD(i,j).*dt(k,:).*dt(l,:)+ ...
    KD(i,l).*dt(k,:).*dt(j,:)+KD(k,j).*dt(i,:).*dt(l,:)+KD(k,l).*dt(i,:).*dt(j,:)) ...
    -60.*dt(i,:).*dt(j,:).*dt(k,:).*dt(l,:)./dr.^7;
                end
            end
        end
    end
    
    % we will contract over our normal vector again

    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:np
            Kijk(i,j,k,l)=normal*squeeze(Psi_ijkl(i,j,k,:,l));
                end
            end
        end
    end
end


function xw = TriGL(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function TriGaussPoints provides the Gaussian points and weights %
% for the Gaussian quadrature of order n for the standard triangles. %
%
%
% Input: n - the order of the Gaussian quadrature (n<=12) %
%
%
% Output: xw - a n by 3 matrix: %
% 1st column gives the x-coordinates of points %
% 2nd column gives the y-coordinates of points %
% 3rd column gives the weights %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xw = zeros(n,3);
if (n == 1)
xw=[0.33333333333333 0.33333333333333 1.00000000000000];
elseif (n == 2)
xw=[0.16666666666667 0.16666666666667 0.33333333333333
0.16666666666667 0.66666666666667 0.33333333333333
0.66666666666667 0.16666666666667 0.33333333333333];
elseif (n == 3)
xw=[0.33333333333333 0.33333333333333 -0.56250000000000
0.20000000000000 0.20000000000000 0.52083333333333
0.20000000000000 0.60000000000000 0.52083333333333
0.60000000000000 0.20000000000000 0.52083333333333];
elseif (n == 4)
xw=[0.44594849091597 0.44594849091597 0.22338158967801
0.44594849091597 0.10810301816807 0.22338158967801
0.10810301816807 0.44594849091597 0.22338158967801
0.09157621350977 0.09157621350977 0.10995174365532
0.09157621350977 0.81684757298046 0.10995174365532
0.81684757298046 0.09157621350977 0.10995174365532];
elseif (n == 5)
xw=[0.33333333333333 0.33333333333333 0.22500000000000
0.47014206410511 0.47014206410511 0.13239415278851
0.47014206410511 0.05971587178977 0.13239415278851
0.05971587178977 0.47014206410511 0.13239415278851
0.10128650732346 0.10128650732346 0.12593918054483
0.10128650732346 0.79742698535309 0.12593918054483
0.79742698535309 0.10128650732346 0.12593918054483];
elseif (n == 6)
xw=[0.24928674517091 0.24928674517091 0.11678627572638
0.24928674517091 0.50142650965818 0.11678627572638
0.50142650965818 0.24928674517091 0.11678627572638
0.06308901449150 0.06308901449150 0.05084490637021
0.06308901449150 0.87382197101700 0.05084490637021
0.87382197101700 0.06308901449150 0.05084490637021
0.31035245103378 0.63650249912140 0.08285107561837
0.63650249912140 0.05314504984482 0.08285107561837
0.05314504984482 0.31035245103378 0.08285107561837
0.63650249912140 0.31035245103378 0.08285107561837
0.31035245103378 0.05314504984482 0.08285107561837
0.05314504984482 0.63650249912140 0.08285107561837];
elseif (n == 7)
xw=[0.33333333333333 0.33333333333333 -0.14957004446768
0.26034596607904 0.26034596607904 0.17561525743321
0.26034596607904 0.47930806784192 0.17561525743321
0.47930806784192 0.26034596607904 0.17561525743321
0.06513010290222 0.06513010290222 0.05334723560884
0.06513010290222 0.86973979419557 0.05334723560884
0.86973979419557 0.06513010290222 0.05334723560884
0.31286549600487 0.63844418856981 0.07711376089026
0.63844418856981 0.04869031542532 0.07711376089026
0.04869031542532 0.31286549600487 0.07711376089026
0.63844418856981 0.31286549600487 0.07711376089026
0.31286549600487 0.04869031542532 0.07711376089026
0.04869031542532 0.63844418856981 0.07711376089026];
elseif (n == 8)
xw=[0.33333333333333 0.33333333333333 0.14431560767779
0.45929258829272 0.45929258829272 0.09509163426728
0.45929258829272 0.08141482341455 0.09509163426728
0.08141482341455 0.45929258829272 0.09509163426728
0.17056930775176 0.17056930775176 0.10321737053472
0.17056930775176 0.65886138449648 0.10321737053472
0.65886138449648 0.17056930775176 0.10321737053472
0.05054722831703 0.05054722831703 0.03245849762320
0.05054722831703 0.89890554336594 0.03245849762320
0.89890554336594 0.05054722831703 0.03245849762320
0.26311282963464 0.72849239295540 0.02723031417443
0.72849239295540 0.00839477740996 0.02723031417443
0.00839477740996 0.26311282963464 0.02723031417443
0.72849239295540 0.26311282963464 0.02723031417443
0.26311282963464 0.00839477740996 0.02723031417443
0.00839477740996 0.72849239295540 0.02723031417443];

else
disp('Bad input n');
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
