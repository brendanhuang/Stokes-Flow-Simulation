% Generate test function f
[X,Y,Z] = meshgrid(1:200,1:100,1:100);
background = sin((X+Y+Z)*pi/100);
blobs = double(sqrt((X-60).^2+(Y-60).^2+(Z-50).^2)<=20) ...
      + double(sqrt((X-140).^2+(Y-60).^2+(Z-50).^2)<=10) ...
      - double(abs(X-100)<=60 & abs(Y-20)<=5 & abs(Z-50)<=20);
f = background + blobs + 0.2*randn(size(X));

% Edge detection parameters
filter = [2; 5];            % Options for Gaussian filter, standard deviation and window size
thresholds = [0.2; 0.1];    % Upper and lower hysteresis thresholds
min_edge_functional = 1e-4; % Low value cutoff for edge detection functional: 
                            % larger value -> faster, less accurate

% Find jump surfaces/edges of f
[edge,conncomp,edgestrength] = edgedetect(f,filter,thresholds,min_edge_functional);

% Segmentation parameters
edgesize = 1;   % Thickness of the voxelized edge surface

% Segment image domain using edge set
seglabels = segment(edge,size(f),edgesize);

% Plot results
figure; imagesc(f(:,:,50)); axis image; set(gca,'YDir','normal');
colormap gray; colorbar; xlabel('x'); ylabel('y'); title('3D-data [z=50]');

figure; title('Edge strength');
plotsurface(edge,edgestrength,0.6); colorbar; view(130,30);

figure; title('Connected components of edge set');
plotsurface(edge,conncomp,0.6); camlight headlight; view(130,30); 

figure; imagesc(seglabels(:,:,50)); axis image; set(gca,'YDir','normal');
xlabel('x'); ylabel('y'); title('Segmentation [z=50]');