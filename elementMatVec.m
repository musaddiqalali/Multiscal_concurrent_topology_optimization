function [keLambda, keMu, feLambda, feMu,KE] = elementMatVec(a, b,Q,select_case)
switch(select_case)
case(1)
DH=Q;
CMu = diag([2 2 1]); CLambda = zeros(3); CLambda(1:2,1:2) = 1;
LM=1;Bele= (1/2/LM)*[-1 0 1 0 1 0 -1 0; 0 -1 0 -1 0 1 0 1;-1 -1 -1 1 1 1 1 -1];
KE=zeros(8,8);
case(2)
DH(1:2,1:2)=Q(1:2,1:2);
CMu = diag([1 1 0]); CLambda = zeros(3);
Bele=0.5*[1,-1,-1,1;1,1,-1,-1];
KE=zeros(4,4);end;
xx=[-1/sqrt(3), 1/sqrt(3)]; yy = xx;
ww=[1,1];
keLambda = zeros(8,8); keMu = zeros(8,8);
feLambda = zeros(8,3); feMu = zeros(8,3);
L = [1 0 0 0; 0 0 0 1; 0 1 1 0];
for ii=1:2
for jj=1:2
x = xx(ii); y = yy(jj);
dNx = 1/4*[-(1-y) (1-y) (1+y) -(1+y)];
dNy = 1/4*[-(1-x) -(1+x) (1+x) (1-x)];
J = [dNx; dNy]*[ -a  a  a  -a;  -b  -b  b  b]';
detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
weight = ww(ii)*ww(jj)*detJ;
G = [invJ zeros(2); zeros(2) invJ];
dN = zeros(4,8);
dN(1,1:2:8) = dNx;
dN(2,1:2:8) = dNy;
dN(3,2:2:8) = dNx;
dN(4,2:2:8) = dNy;
B = L*G*dN;
KE = KE + ww(ii)*ww(jj)*detJ*Bele'*DH*Bele;
keLambda = keLambda + weight*(B' * CLambda * B);
keMu = keMu + weight*(B' * CMu * B);
feLambda = feLambda + weight*(B' * CLambda * diag([1 1 1]));
feMu = feMu + weight*(B' * CMu * diag([1 1 1]));
end;end;end