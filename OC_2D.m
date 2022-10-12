function [x, xPhys, change] = OC_2D(x, dc,volfrac,rmin, nelx,nely)
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));sH = zeros(size(iH));k = 0;
for i1 = 1:nelx
for j1 = 1:nely
e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);dv = ones(nely,nelx);
%% Filtering
dc(:) = H*(dc(:)./Hs);
dv(:) = H*(dv(:)./Hs);
%% Optimality Criteria Update of Design Variables & Physical Densities
l1 = 0; l2 = 1e9; move = 0.009;
while (l2-l1)/(l1+l2) > 1e-6
lmid = 0.5*(l2+l1);
xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(abs(dc)./dv/lmid)))));
xPhys(:) = (H*xnew(:))./Hs;
if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
end;  change = max(abs(xnew(:)-x(:)));x = xnew;end