%% Concurrent Multiphysics Multiscale Topology Optimization 
%% Material input
%% Microscale input
Micro_nelx=100;Micro_nely=100;Micro_volfrac=.5;penal=3;Micro_rmin=1.5;
% Macroscale input
Macro_nelx=200;Macro_nely=90;Macro_volfrac=.5;%Secound Example  Macro_nelx=150;Macro_nely=150
Macro_rmin=1.5;dx = 1/Macro_nelx; dy = 1/Macro_nely;
% Optimization input
maxIter=500;eta=0.5;
%% Mechanical FEM
Macro_nodenrs_mech = reshape(1:(1+Macro_nelx)*(1+Macro_nely),1+Macro_nely,1+Macro_nelx);
Macro_edofVec_mech = reshape(2*Macro_nodenrs_mech(1:end-1,1:end-1)+1,Macro_nelx*Macro_nely,1);
Macro_edofMat_mech = repmat(Macro_edofVec_mech,1,8)+repmat([0 1 2*Macro_nely+[2 3 0 1] -2 -1], ...
Macro_nelx*Macro_nely,1);
iK_mech = reshape(kron(Macro_edofMat_mech,ones(8,1))',64*Macro_nelx*Macro_nely,1);
jK_mech = reshape(kron(Macro_edofMat_mech,ones(1,8))',64*Macro_nelx*Macro_nely,1);
U = zeros(2*(Macro_nely+1)*(Macro_nelx+1),1);
F = zeros(2*(Macro_nely+1)*(Macro_nelx+1),1);
F(2*(Macro_nelx+1)*(Macro_nely+1),1)=1 ;  
fixeddofs_mech=union([2:2*(Macro_nely+1)],[1]);
alldofs_mech = [1:2*(Macro_nely+1)*(Macro_nelx+1)];
freedofs_mech = setdiff(alldofs_mech,fixeddofs_mech);
%% Heat FEM
Macro_nodenrs = reshape(1:(1+Macro_nelx)*(1+Macro_nely),1+Macro_nely,1+Macro_nelx);
Macro_edofVec = reshape(Macro_nodenrs(1:end-1,1:end-1)+1,Macro_nelx*Macro_nely,1);
Macro_edofMat_Heat = repmat(Macro_edofVec,1,4)+repmat([0 Macro_nely+[1 0] -1],Macro_nelx* ...
Macro_nely,1);
iK = reshape(kron(Macro_edofMat_Heat,ones(4,1))',16*Macro_nelx*Macro_nely,1);
jK = reshape(kron(Macro_edofMat_Heat,ones(1,4))',16*Macro_nelx*Macro_nely,1);
P = sparse((Macro_nely+1)*(Macro_nelx+1),1); P(:,1) = .01;
fixeddofs_Heat = [Macro_nely/2-20:Macro_nely/2+20];
T = zeros((Macro_nely+1)*(Macro_nelx+1),1);
alldofs_Heat = [1:(Macro_nely+1)*(Macro_nelx+1)];
freedofs_Heat = setdiff(alldofs_Heat,fixeddofs_Heat);
%% Mechanical  Compliance **
[Macro_x,Micro_x]=get_Initials(Macro_volfrac,Micro_volfrac,Micro_nelx,Micro_nely, ...
Macro_nelx,Macro_nely);
Micro_xPhys = Micro_x;Macro_xPhys=Macro_x;
[QM, ~] = Homogenization_full(Micro_xPhys,1);
[~,~,~,~,KE_Homogenized_Mech] = elementMatVec(dx/2, dy/2,QM,1);
[N_macro_Mech]=Objective_Calculator(iK_mech,jK_mech,Macro_nely,Macro_nelx,freedofs_mech, ...
Macro_edofMat_mech,Macro_xPhys,penal,KE_Homogenized_Mech,F,1);
[U_macro_Mech,~]=Concurrent_solver(Micro_nelx,Micro_nely,Macro_nelx,Macro_nely, ...
Macro_edofMat_mech,iK_mech,jK_mech,F,U,freedofs_mech,freedofs_Heat,Macro_edofMat_Heat, ...
T,P,iK,jK,Macro_volfrac,Micro_volfrac,Micro_rmin,Macro_rmin,dx,dy,penal,maxIter,1);
Diffrence_Nadir_Utopia_Mech=N_macro_Mech-U_macro_Mech;
clear Macro_x Micro_x
%% Heat Compliance **
[Macro_x,Micro_x]=get_Initials(Macro_volfrac,Micro_volfrac,Micro_nelx,Micro_nely,Macro_nelx,Macro_nely);
Micro_xPhys = Micro_x;Macro_xPhys=Macro_x;
[Q,~]=Homogenization_full(Micro_xPhys,2);
[~,~,~,~,KE_Homogenized_Heat] = elementMatVec(dx/2, dy/2,Q,2);
[N_macro_Heat]=Objective_Calculator(iK,jK,Macro_nely,Macro_nelx,freedofs_Heat,Macro_edofMat_Heat, ...
Macro_xPhys,penal,KE_Homogenized_Heat,P,2);
[U_macro_Heat,~]=Concurrent_solver(Micro_nelx,Micro_nely,Macro_nelx,Macro_nely,Macro_edofMat_mech, ...
iK_mech,jK_mech,F,U,freedofs_mech,freedofs_Heat,Macro_edofMat_Heat,T,P,iK,jK,Macro_volfrac, ...
Micro_volfrac,Micro_rmin,Macro_rmin,dx,dy,penal,maxIter,2);
Diffrence_Nadir_Utopia_Heat=N_macro_Heat-U_macro_Heat;
%% MOO optimization loop
clear Macro_x Micro_x
[Macro_x,Micro_x]=get_Initials(Macro_volfrac,Micro_volfrac,Micro_nelx,Micro_nely,Macro_nelx,Macro_nely);
Micro_xPhys = Micro_x;Macro_xPhys=Macro_x;loop = 0;
while loop < maxIter;
loop = loop+1; 
%% Mechanical Part
[QM,dQM]=Homogenization_full(Micro_xPhys,1);
[~,~,~,~,KE_Macro_Homogenized_mech] = elementMatVec(dx/2, dy/2,QM,1);
sK_mech = reshape(KE_Macro_Homogenized_mech(:)*(1e-6+Macro_xPhys(:)'.^penal*(1-1e-6)),64*Macro_nelx*Macro_nely,1);
K_mech = sparse(iK_mech,jK_mech,sK_mech); K_mech = (K_mech+K_mech')/2;
U(freedofs_mech,:) = K_mech(freedofs_mech,freedofs_mech)\F(freedofs_mech,:);
ce_mech = reshape(sum((U(Macro_edofMat_mech)*KE_Macro_Homogenized_mech).*U(Macro_edofMat_mech),2),Macro_nely,Macro_nelx);
c= sum(sum((1e-6+Macro_xPhys.^penal*(1-1e-6)).*ce_mech));
Macro_dc_mech = -penal*(1-1e-6)*Macro_xPhys.^(penal-1).*ce_mech;
Micro_dc_mech = zeros(Micro_nely, Micro_nelx);
for i = 1:Micro_nelx*Micro_nely
dQe_mech = [dQM{1,1}(i) dQM{1,2}(i) dQM{1,3}(i);
            dQM{2,1}(i) dQM{2,2}(i) dQM{2,3}(i);
            dQM{3,1}(i) dQM{3,2}(i) dQM{3,3}(i)];
[~,~,~,~,dKE_mech] = elementMatVec(dx/2, dy/2,dQe_mech,1);
dce_mech = reshape(sum((U(Macro_edofMat_mech)*dKE_mech).*U(Macro_edofMat_mech),2),Macro_nely,Macro_nelx);
Micro_dc_mech(i) = -sum(sum((1e-6+Macro_xPhys.^penal*(1-1e-6)).*dce_mech));
end
%% Heat Part
[Q,dQ]=Homogenization_full(Micro_xPhys,2);
[~,~,~,~,KE_Homogenized_Heat] = elementMatVec(dx/2, dy/2,Q,2);
T = zeros((Macro_nely+1)*(Macro_nelx+1),1);
sK = reshape(KE_Homogenized_Heat(:)*(1e-6+Macro_xPhys(:)'.^penal*(1-1e-6)),16*Macro_nelx*Macro_nely,1);
KH = sparse(iK,jK,sK); KH = (KH+KH')/2;
T(freedofs_Heat,:) = KH(freedofs_Heat,freedofs_Heat)\P(freedofs_Heat,:);
ce_Heat = reshape(sum((T(Macro_edofMat_Heat)*KE_Homogenized_Heat).*T(Macro_edofMat_Heat),2),Macro_nely,Macro_nelx);
c_Heat = sum(sum((1e-6+Macro_xPhys.^penal*(1-1e-6)).*ce_Heat)) 
Macro_dcH = -penal*(1-1e-6)*Macro_xPhys.^(penal-1).*ce_Heat;
Macro_dv = ones(Macro_nely, Macro_nelx);
Micro_dcH = zeros(Micro_nely, Micro_nelx);
for i = 1:Micro_nelx*Micro_nely;  
dQe = [dQ{1,1}(i) dQ{1,2}(i);
dQ{2,1}(i) dQ{2,2}(i)];
[~,~,~,~,dKE] = elementMatVec(dx/2, dy/2,dQe,2);
dce = reshape(sum((T(Macro_edofMat_Heat)*dKE).*T(Macro_edofMat_Heat),2),Macro_nely,Macro_nelx);
Micro_dcH(i) = -sum(sum((1e-6+Macro_xPhys.^penal*(1-1e-6)).*dce));
end
Micro_dv = ones(Micro_nely, Micro_nelx);
%% MMO objective function
Macro_dc= eta*Macro_dcH/(Diffrence_Nadir_Utopia_Heat)+(1-eta)*Macro_dc_mech/(Diffrence_Nadir_Utopia_Mech);
Micro_dc=eta*Micro_dcH/(Diffrence_Nadir_Utopia_Heat)+(1-eta)*Micro_dc_mech/(Diffrence_Nadir_Utopia_Mech);
%% Optimality Criteria Update for Macro and Micro Element Densities
[Macro_x, Macro_xPhys, Macro_change] = OC_2D(Macro_x, Macro_dc, Macro_volfrac,Macro_rmin, Macro_nelx,Macro_nely);
[Micro_x, Micro_xPhys, Micro_change] = OC_2D(Micro_x, Micro_dc, Micro_volfrac,Micro_rmin, Micro_nelx,Micro_nely);
Macro_xPhys = reshape(Macro_xPhys, Macro_nely, Macro_nelx); Micro_xPhys = reshape(Micro_xPhys, Micro_nely, Micro_nelx);
%% Printing the Results
fprintf(' It.:%5i Obj.:%11.4f Macro_Vol.:%7.3f Micro_Vol.:%7.3f Macro_ch.:%7.3f Micro_ch.:%7.3f4\n',loop,mean(Macro_xPhys(:)),...
mean(Micro_xPhys(:)), Macro_change, Micro_change);
hold on;
figure (6);clf;colormap(gray); imagesc(1-Micro_x); caxis([0 1]); 
axis equal; axis off; t1=title('Microscale Design');t1.Color ='#FE8402';drawnow;
figure (7);clf;colormap(gray); imagesc(1-Macro_x); caxis([0 1]);
axis equal; axis off; t2=title('Macroscale Design');t2.Color ='#0053ED';drawnow;
end;
disp('********* Concurrent Multiphysics Multiscale Topology Optimization is finished ***************');
