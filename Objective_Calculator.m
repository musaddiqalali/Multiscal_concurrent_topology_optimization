function [c]=Objective_Calculator(iK,jK,nely,nelx,freedofs,edofMat,x,penal,KE,F,option)
switch(option)
case(1) %Mechanical   
U = zeros(2*(nely+1)*(nelx+1),1);
sK = reshape(KE(:)*(1e-6+x(:)'.^penal*(1-1e-6)),64*nely*nelx,1);K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
ce= reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
c= sum(sum((1e-6+x.^penal*(1-1e-6)).*ce ));  
case(2) %Thermal
U = zeros((nely+1)*(nelx+1),1);
sK = reshape(KE(:)*(1e-6+x(:)'.^penal*(1-1e-6)),16*nelx*nely,1);K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
c = sum(sum((1e-6+x.^penal*(1-1e-6)).*ce)); 
end