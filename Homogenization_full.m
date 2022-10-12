function [CH,dCH]= Homogenization_full(Micro,select_case)
%% Initialize input data
[nelx,nely]=size(Micro);
for i=1:nelx
for j=1:nely
if Micro(i,j)<.1;x(i,j)=2;else;x(i,j)=1;
end;end;end
lx=1;ly=1;lambda=[1 0];mu=[1 0];
dx = lx/nelx; dy = ly/nely;
nel = nelx*nely;
Q=ones(3,3);
[keLambda, keMu, feLambda, feMu,~] = elementMatVec(dx/2, dy/2,Q,select_case);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nel,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nel,1);
switch(select_case);case(1);
case(2) 
keMu(1:2:end,1:2:end) = keMu(1:2:end,1:2:end)+keMu(2:2:end, 2:2:end);end
%% Impose Periodic Boundary Conditions
nn = (nelx+1)*(nely+1); 
nnP = (nelx)*(nely); 
nnPArray = reshape(1:nnP, nely, nelx);
nnPArray(end+1,:) = nnPArray(1,:);
nnPArray(:,end+1) = nnPArray(:,1);
dofVector = zeros(2*nn, 1);
dofVector(1:2:end) = 2*nnPArray(:)-1;
dofVector(2:2:end) = 2*nnPArray(:);
edofMat = dofVector(edofMat);
ndof = 2*nnP;
%% Assemble Stiffness Matrix
iK = kron(edofMat,ones(8,1))';
jK = kron(edofMat,ones(1,8))';
lambda = lambda(1)*(x==1) + lambda(2)*(x==2);
mu = mu(1)*(x==1) + mu(2)*(x==2);
for i=1 nely;for j=1;
lambda(i,j)=lambda(i,j);
mu(i,j)=mu(i,j);end;end
sK = keLambda(:)*lambda(:).' + keMu(:)*mu(:).';
K = sparse(iK(:), jK(:), sK(:), ndof, ndof);
%% Load Vectors and Solution
sF = feLambda(:)*lambda(:).'+feMu(:)*mu(:).';
iF = repmat(edofMat',3,1);
jF = [ones(8,nel); 2*ones(8,nel); 3*ones(8,nel)];
F = sparse(iF(:), jF(:), sF(:), ndof, 3);
% Solve 
activedofs=edofMat(x==1,:);
activedofs=sort(unique(activedofs(:)));
switch(select_case);
case(1) 
Xi=zeros(ndof,3);
Xi(activedofs(3:end),:) =K(activedofs(3:end),activedofs(3:end))\F(activedofs(3:end),:);
case(2)
Xi = zeros(ndof,2);
Xi(activedofs(3:2:end),:) = K(activedofs(3:2:end),activedofs(3:2:end))\...
[F(activedofs(3:2:end),1) F(activedofs(4:2:end),2)];end
%% Homogenization
Xi0 = zeros(nel, 8, 3);Xi0_e = zeros(8, 3);
ke = keMu + keLambda;fe = feMu + feLambda;
Loop=0;
switch(select_case);
case(1)
Loop=3;
Xi0_e([3 5:end],:) = ke([3 5:end],[3 5:end])\fe([3 5:end],:);
case(2)
Loop=2;
Xi0_e(3:2:end,1:2) = keMu(3:2:end,3:2:end)\[feMu(3:2:end,1) feMu(4:2:end,2)]; end
Xi0(:,:,1) = kron(Xi0_e(:,1)', ones(nel,1));
Xi0(:,:,2) = kron(Xi0_e(:,2)', ones(nel,1));
Xi0(:,:,3) = kron(Xi0_e(:,3)', ones(nel,1));
CH = zeros(3);
cellVolume = lx*ly;
for i = 1:Loop
for j = 1:Loop
sumLambda = ((Xi0(:,:,i) - Xi(edofMat+(i-1)*ndof))*keLambda).*...
(Xi0(:,:,j) - Xi(edofMat+(j-1)*ndof));
sumMu = ((Xi0(:,:,i) - Xi(edofMat+(i-1)*ndof))*keMu).*...
(Xi0(:,:,j) - Xi(edofMat+(j-1)*ndof));
sumLambda = reshape(sum(sumLambda,2), nely, nelx);
sumMu = reshape(sum(sumMu,2), nely, nelx);
qe=lambda.*sumLambda + mu.*sumMu;
CH(i,j) =  1/cellVolume*sum(sum(qe.*(1e-6+Micro.^(3)*(1-1e-6))));
dCHq = 1/cellVolume*(qe.*(3*Micro.^(3-1)*(1-0.0001)).'); dCH{i,j}=dCHq;
end;end