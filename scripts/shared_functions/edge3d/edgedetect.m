function [edge,conn_comp,edge_strength] = edgedetect(f,filter,thresholds,min_edge_functional)
% function [edge,conn_comp,edge_strength] = edgedetect(f,filter,thresholds,min_edge_functional)
%
% Detect edges/jumps of given 3D-image f using Canny edge detection
% (differential version, cf. Lindeberg 1998) with hysteresis thresholding
%
% Output: 
%   edge: Edge/jump surfaces (face & vertex list)
%   conn_comp: Connected components of edge set
%   edge_strength: Gradient magnitude |nabla f| (of smoothed data)
%
% Input:
%   f: Data to detect jumps in (3D-array containing grayscale values)
%   filter: (2x1)-vector, filter(1) contains the standard deviation and
%           filter(2) the window size of Gaussian convolution filter
%   thresholds: (2x1)-vector, upper/lower hysteresis thresholds.
%               Edge parts with gradient magnitude less then thresholds(2)
%               are removed from the set set. Only those connected
%               componenents where the gradient magnitude rises above
%               thresholds(1) are kept in the edge set.
%   min_edge_functional: Cutoff for edge detection differential for
%                        speedup: larger value -> faster, less accurate.
%                        Set to 0 if no speedup is required. 

% default input
if nargin <= 3
  min_edge_functional=0;
end
if nargin <= 2
  thresholds=0;
end
if nargin <= 1
  filter=0;  
end
if nargin == 0
  error('Missing input arguments.');
end

% verify input arguments
if ndims(f)~=3 || ~isvector(thresholds) || ~isvector(filter) ...
               || ~isscalar(min_edge_functional)
    error('Wrong input argument structure.');
end

% standard window size: 2*sigma
if isscalar(filter)
    window=ceil(2*filter);
else
    window=filter(2); sigma=filter(1);
end

% apply smoothing filter
if(sigma > 0)    
    filt=gaussfilter([window;window;window],(diag(sigma).^2));
    f_smooth=imfilter(f,filt,NaN);
else
    f_smooth=f;
end

% calculate gradient magnitude
[gx,gy,gz]=gradient(f_smooth);
f_gradmag=sqrt(gx.^2+gy.^2+gz.^2);
    
% define differential filter masks - first order
filt_dx=[-1/2, 0, 1/2];
filt_dy=shiftdim(filt_dx,1);
filt_dz=shiftdim(filt_dx,-1);

% define differential filter masks - second order
filt_dxx=[1, -2, 1];
filt_dyy=shiftdim(filt_dxx,1);
filt_dzz=shiftdim(filt_dxx,-1);
filt_dxy=[1/4, 0, -1/4; 0 0 0; -1/4, 0, 1/4];
filt_dyz=shiftdim(filt_dxy,-1);
filt_dxz=shiftdim(filt_dyz,1);

% define differential filter masks - third order
filt_dxxx=[-1/2, 1, 0, -1, 1/2];
filt_dyyy=shiftdim(filt_dxxx,1);
filt_dzzz=shiftdim(filt_dxxx,-1);
filt_dxxy=[-1/2, 1, -1/2; 0 0 0; 1/2, -1, 1/2];
filt_dxxz=permute(filt_dxxy,[3 2 1]);
filt_dyyx=permute(filt_dxxy,[2 1 3]);
filt_dyyz=permute(filt_dxxy,[2 3 1]);
filt_dzzx=permute(filt_dxxy,[3 1 2]);
filt_dzzy=permute(filt_dxxy,[1 3 2]);
filt_dxyz=(1/2)*cat(3,-filt_dxy,zeros(3,3),filt_dxy);

% apply filters
f_dx=imfilter(f_smooth,filt_dx);
f_dy=imfilter(f_smooth,filt_dy);
f_dz=imfilter(f_smooth,filt_dz);
f_dxx=imfilter(f_smooth,filt_dxx);
f_dyy=imfilter(f_smooth,filt_dyy);
f_dzz=imfilter(f_smooth,filt_dzz);
f_dxy=imfilter(f_smooth,filt_dxy);
f_dyz=imfilter(f_smooth,filt_dyz);
f_dxz=imfilter(f_smooth,filt_dxz);
f_dxxx=imfilter(f_smooth,filt_dxxx);
f_dyyy=imfilter(f_smooth,filt_dyyy);
f_dzzz=imfilter(f_smooth,filt_dzzz);
f_dxxy=imfilter(f_smooth,filt_dxxy);
f_dxxz=imfilter(f_smooth,filt_dxxz);
f_dyyx=imfilter(f_smooth,filt_dyyx);
f_dyyz=imfilter(f_smooth,filt_dyyz);
f_dzzx=imfilter(f_smooth,filt_dzzx);
f_dzzy=imfilter(f_smooth,filt_dzzy);
f_dxyz=imfilter(f_smooth,filt_dxyz);

% calculate differential used for non-maximum surpression
maxedge = f_dx.^2.*f_dxx + f_dy.^2.*f_dyy  ...
    + f_dz.^2.*f_dzz + 2*f_dx.*f_dy.*f_dxy ...
    + 2*f_dy.*f_dz.*f_dyz + 2*f_dx.*f_dz.*f_dxz;

% remove small variation of edge detection functional
maxedge(abs(maxedge)<min_edge_functional)=0;

% find potential edges by taking zero pass of above function
edge=isosurface(maxedge,0);  

% return if threshold is zero or edge set is empty
if ~thresholds(1) || isempty(edge.vertices)
    return; 
end

% third order sign condition to filter out extrema which are not maxima
signcond = f_dx.^3.*f_dxxx + f_dy.^3.*f_dyyy + f_dz.^3.*f_dzzz ...
         + 3*f_dx.^2.*f_dy.*f_dxxy + 3*f_dx.^2.*f_dz.*f_dxxz ...
         + 3*f_dy.^2.*f_dx.*f_dyyx + 3*f_dy.^2.*f_dz.*f_dyyz ...
         + 3*f_dz.^2.*f_dx.*f_dzzx + 3*f_dz.^2.*f_dy.*f_dzzy ...
         + 6*f_dx.*f_dy.*f_dz.*f_dxyz;

% apply third order sign condition 
edge_sign=interp3(signcond,edge.vertices(:,1),edge.vertices(:,2),edge.vertices(:,3),'linear');
ind=(edge_sign<0); edge=reducesurface(edge,ind);

% interpolate gradient magnitude on edge vertices for thresholding
edge_strength=interp3(f_gradmag,edge.vertices(:,1),edge.vertices(:,2),edge.vertices(:,3),'linear');

% remove edge vertices with gradient magnitude below lower threshold
ind=(edge_strength>thresholds(2)); edge=reducesurface(edge,ind); edge_strength=edge_strength(ind);

% remove connected components of edge set where gradient magnitude 
% is never above the upper threshold
[ind,conn_comp]=hysteresis(edge,edge_strength,thresholds);
edge=reducesurface(edge,ind); edge_strength=edge_strength(ind); 
[~,~,conn_comp]=unique(conn_comp(ind));

end


function [ind,conn_comp] = hysteresis(edge,edge_strength,thresholds)
% Apply upper hysteresis threshold on face & vertex list 

% initialize accepted vertex list
ind=true(size(edge.vertices,1),1);

% find connected components of edge set
conn_comp=conncomp(edge);

% loop over connected components
for i=1:max(conn_comp)
    
    % remove component if maximal edge strength is below upper threshold
    if max(edge_strength(conn_comp==i)) < thresholds(1)
        ind(conn_comp==i)=false;
    end
    
end

end

function surf_new = reducesurface(surf,ind)
% Reduce list of vertices and faces (triangles) to given indices

% reduce vertex list
surf_new.vertices=surf.vertices(ind,:);

% corresponding faces
num=(1:sum(ind))'; indlist=NaN(size(surf.vertices,1),1);
indlist(ind)=num; new_elem=indlist(surf.faces);
surf_new.faces=new_elem(~isnan(sum(new_elem,2)),:);

end

function cc = conncomp(surf,label)
% Find connected components of given surface (list of faces and vertices)

% initialize label list if not given a prior
if ~exist('label','var')
   label=(1:size(surf.vertices,1))';
end

% loop until all nodes of every face have the same label
while max(std(label(surf.faces),0,2))~=0
    
    % smallest label of every face (vectorized)
    m=repmat(min(label(surf.faces),[],2),[3 1]);

    % set label of node to smallest of face (randomize to break cycles)
    j=surf.faces(:); p=randperm(length(m));
    label(j(p))=m(p);    
end

% relabel
[~,~,cc]=unique(label);

end

function filter = gaussfilter(sz,sigma)
% Calculate Gaussian filter tensor with given dimensions sz
% and covariance matrix sigma

% initialize
dims=length(sz); ind_cell=cell(dims,1); filter=zeros(sz');

% calculate midpoint
midpoint=(sz+1)/2;

% define (unnormalized) gaussian function
gauss=@(x) exp(-x'*(sigma\x)/2);

% loop over all elements of tensor
for i=1:numel(filter)
    
    % get indices of element with linear index i
    [ind_cell{:}]=ind2sub(sz',i); ind=cell2mat(ind_cell);
    
    % calculate gaussian value at size i
    filter(i)=gauss(ind-midpoint);
    
end

% get l1 norm of filter
n=filter;

for d=1:dims
    n=sum(n,dims-d+1);
end

% normalize
filter=filter/n;

end