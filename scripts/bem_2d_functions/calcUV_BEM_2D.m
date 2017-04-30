function [u,v,pres,stress] = calcUV_BEM_2D(f0,u0,SLmat,DLmat,xpts,ypts,visc,PSL,PDL,SSL,SDL)
    
uvtot = -1/(4*pi*visc)*SLmat*f0+1/(4*pi)*DLmat*u0;


[xgrid,ygrid]=ndgrid(linspace(1,xpts,xpts),linspace(1,ypts,ypts));
u1=uvtot(1:xpts*ypts);
v1=uvtot(xpts*ypts+1:end);

u=accumarray([xgrid(:),ygrid(:)],u1,[xpts ypts]);
v=accumarray([xgrid(:),ygrid(:)],v1,[xpts ypts]);

if nargout==3
    pres1 = -1/(4*pi)*PSL*f0+visc/(4*pi)*PDL*u0;
    pres=accumarray([xgrid(:),ygrid(:)],pres1,[xpts ypts]); 
elseif nargout ==4
    pres1 = -1/(4*pi)*PSL*f0+visc/(4*pi)*PDL*u0;
    stress1 = -1/(4*pi)*SSL*f0+visc/(4*pi)*SDL*u0;
    pres=accumarray([xgrid(:),ygrid(:)],pres1,[xpts ypts]); 
    stress1=reshape(stress1,[xpts*ypts,2,2]);
    stress=zeros(xpts,ypts,2,2);
    for i=1:2
       for j=1:2
           stress(:,:,i,j)=accumarray([xgrid(:),ygrid(:)],squeeze(stress1(:,i,j)),[xpts ypts]);
       end
    end
end

end