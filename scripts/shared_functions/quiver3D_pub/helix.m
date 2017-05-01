function [pos,mag] = helix(radius, height, numRevolutions, numPoints, arrowScale)

% [pos,mag] = helix(radius, height, numRevolutions, numPoints, arrowScale) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     Parametric equation of the 3D helix
%         [pos] = helix(radius, height) - output pos = [x(:),y(:),z(:)] locations along a helix
%         [pos,mag] = helix(radius, height) - output (x,y,z) along with tangent arrows (u,v,w) with fixed magnitude
%         [...] = helix(..., numRevolutions) - default = two rotations of the helix
%         [...] = helix(..., numRevolutions, numPoints) - number of samples along the helix (default = 25)
%         [...] = helix(..., numRevolutions, numPoints, arrowScale) - magnitude of the tangent arrows (default = 0.8)
%     
% Example:
%    [pos, mag] = helix(5, 4, 7, 50, 0.24);  
%    quiver3D(posArray1, magnitudeArray1);   %---- plots a helix with radius=5, height=4, with 7 revolutions, sampled at 50 points with arrows 24% of the length
%                                            %     between sampled points
% 
%     Author: Shawn Arseneau
%     Created: September 15, 2006
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin<2 || nargin>5         
        error('Incorrect number of inputs to helix');   
    end
    if nargout~=1 && nargout~=2     
        error('Incorrect number of outputs to helix');  
    end
    if nargin<=4                    
        arrowScale = 0.8;                               
    end
    if nargin<=3                    
        numPoints = 25;                                 
    end
    if nargin<=2                    
        numRevolutions = 2;                               
    end

    endAngle = numRevolutions*2*pi;
    t = (0:endAngle/numPoints:endAngle);
    X_all = radius*cos(t);
    Y_all = radius*sin(t);
    Z_all = height*t;
    
    X = X_all(1:end-1);   X_end = X_all(2:end);
    Y = Y_all(1:end-1);   Y_end = Y_all(2:end);
    Z = Z_all(1:end-1);   Z_end = Z_all(2:end);
    
    pos = [X(:),Y(:),Z(:)];
    
    if nargout==2
        U = arrowScale*(X_end - X);
        V = arrowScale*(Y_end - Y);
        W = arrowScale*(Z_end - Z);
        mag = [U(:),V(:),W(:)];
    end



