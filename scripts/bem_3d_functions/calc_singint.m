%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stokes-Flow-Simulation by Brendan Huang
% Script: calc_singint
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
% Description: This function analytically calculates the Green's function
% (single layer potential) between a right triangular boundary element and
% a point that is located at the center of mass of that triangle. The
% function breaks up the integration into three segments corresponding to 
% three vertices and evaluates the integral. Specifically, it uses a radial
% transform of coordinates that takes the location of the singularity to be
% at r=0
% 
%
% Inputs:
% r1, r2, r3 - vertices of the triangle
% p - center of mass of the triangle 
%
% Outputs: 
% Gxx ... Gzz: Nine components of the single layer potential 
%
% References: I derived the analytic expression for each of the three
% segments of the triangle.  I could not find any references that published
% the form here

function [Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,Gzx,Gzy,Gzz] = calc_singint(p,r1,r2,r3)

% First we need to calculate a few properties of the triangle, including the
% lengths of the sides, the vector between the center of mass and each
% point
s12=norm(r1-r2);
s13=norm(r1-r3);
s23=norm(r2-r3);
Avec=r1-p;
Bvec=r2-p;
Cvec=r3-p;
A=norm(Avec);
B=norm(Bvec);
C=norm(Cvec);

% check clockwise versus counterclockwise orientation of triangle. If it is
% oriented clockwise, we need to reverse direction of integration.
c1=cross([r2-r1;0],[r3-r1;0]);

% get our angular limits of integration
theta01=atan2(Avec(2),Avec(1));
theta02=atan2(Bvec(2),Bvec(1));
theta03=atan2(Cvec(2),Cvec(1));
arc1=theta02-theta01;
    if abs(arc1)>pi
        arc1=abs(abs(arc1)-2*pi);
    else
        arc1=abs(arc1);
    end
    
arc2=theta03-theta02;

    if abs(arc2)>pi
        arc2=abs(abs(arc2)-2*pi);
    else
        arc2=abs(arc2);
    end
    
arc3=theta01-theta03;

    if abs(arc3)>pi
        arc3=abs(abs(arc3)-2*pi);
    else
        arc3=abs(arc3);
    end

% use law of cosines to calculate angles of triangle
thetaA1=acos((B^2+s12^2-A^2)/(2*B*s12));
thetaB1=acos((A^2+s12^2-B^2)/(2*A*s12));
thetaB2=acos((C^2+s23^2-B^2)/(2*C*s23));
thetaC2=acos((B^2+s23^2-C^2)/(2*B*s23));
thetaC3=acos((A^2+s13^2-C^2)/(2*A*s13));
thetaA3=acos((C^2+s13^2-A^2)/(2*C*s13));

%set limits of integration
if c1(3)>0
    theta1s=theta01;
    theta2s=theta02;
    theta3s=theta03;
    theta1a=thetaB1;
    theta2a=thetaC2;
    theta3a=thetaA3;
    AA=A;
    BB=B;
    CC=C;
else
    theta1s=theta02;
    theta2s=theta03;
    theta3s=theta01;
    theta1a=thetaA1;
    theta2a=thetaB2;
    theta3a=thetaC3;
    AA=B;
    BB=C;
    CC=A;
end

%evaluate integrals analytically we have to compute each of the vectorial
%components (coupling xx, xy, xz etc.). Amongst each component, each
%integral consists of three segments, arising from the three sides of the
%triangle
funintxy= @(phi,A,b,phi0) 1/2*A*sin(b)*((log(cos(1/2*(b-phi0+phi))) ...
    -log(sin(1/2*(b-phi0+phi))))*sin(2*(b-phi0))-2*sin(b-phi0-phi));
intxy=feval(funintxy,theta1s+arc1,AA,theta1a,theta1s)- feval(funintxy,theta1s,AA,theta1a,theta1s) ...
    + feval(funintxy,theta2s+arc2,BB,theta2a,theta2s)- feval(funintxy,theta2s,BB,theta2a,theta2s) ...
    + feval(funintxy,theta3s+arc3,CC,theta3a,theta3s)- feval(funintxy,theta3s,CC,theta3a,theta3s);

funintxx= @(phi,A,b,phi0) A*sin(b)*(cos(b - phi0)*cos(phi) + ...
   cos(b - phi0).^2*(-log(cos(1/2*(b - phi0 + phi))) + ...
      log(sin(1/2*(b - phi0 + phi)))) + sin(b - phi0)*sin(phi));
intxx=feval(funintxx,theta1s+arc1,AA,theta1a,theta1s)- feval(funintxx,theta1s,AA,theta1a,theta1s) ...
    + feval(funintxx,theta2s+arc2,BB,theta2a,theta2s)- feval(funintxx,theta2s,BB,theta2a,theta2s) ...
    + feval(funintxx,theta3s+arc3,CC,theta3a,theta3s)- feval(funintxx,theta3s,CC,theta3a,theta3s);

funintyy= @(phi,A,b,phi0) 1/2*A*sin(b)*(-2*cos(b - phi - phi0) + ... 
      2*(-log(cos(1/2*(b + phi - phi0))) + ... 
      log(sin(1/2*(b + phi - phi0))))*sin(b - phi0)^2);  
intyy=feval(funintyy,theta1s+arc1,AA,theta1a,theta1s)- feval(funintyy,theta1s,AA,theta1a,theta1s) ...
    + feval(funintyy,theta2s+arc2,BB,theta2a,theta2s)- feval(funintyy,theta2s,BB,theta2a,theta2s) ...
    + feval(funintyy,theta3s+arc3,CC,theta3a,theta3s)- feval(funintyy,theta3s,CC,theta3a,theta3s);

funintr= @(phi,A,b,phi0) A*(-log(cos(1/2*(b + phi - phi0))) ... 
   + log(sin(1/2*(b + phi - phi0))))*sin(b);  
intr=feval(funintr,theta1s+arc1,AA,theta1a,theta1s)- feval(funintr,theta1s,AA,theta1a,theta1s) ...
    + feval(funintr,theta2s+arc2,BB,theta2a,theta2s)- feval(funintr,theta2s,BB,theta2a,theta2s) ...
    + feval(funintr,theta3s+arc3,CC,theta3a,theta3s)- feval(funintr,theta3s,CC,theta3a,theta3s);

intzz=0;
intyz=0;
intxz=0;

% Get matrix coefficients based on Gij = dij/r + xi xj/r^3
Gxx=intxx+intr;
Gxy=intxy;
Gxz=intxz;
Gyx=intxy;
Gyy=intyy+intr;
Gyz=intyz;
Gzx=intxz;
Gzy=intxz;
Gzz=intzz+intr;


end