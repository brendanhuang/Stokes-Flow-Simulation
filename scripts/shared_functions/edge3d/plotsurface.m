function plotsurface(surf,coloring,alpha)
% function plotsurface(surf,coloring,alpha)
% Plot surface (face & vertex list) with given coloring

% default: uniform coloring
if ~exist('coloring','var')
    coloring=ones(length(surf.vertices),1);
end

% default: facealpha = 0.5
if ~exist('alpha','var');
    alpha=0.5;
end

patch('Faces',surf.faces,'Vertices',surf.vertices,'FaceVertexCData',coloring,'FaceColor','interp','edgecolor', 'interp','FaceAlpha',alpha); 
xlabel('x'); ylabel('y'); zlabel('z'); view(3); daspect([1 1 1]);
shading flat;

end

