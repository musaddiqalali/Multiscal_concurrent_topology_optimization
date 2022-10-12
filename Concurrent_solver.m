function [C_concurrent,CH]=Concurrent_solver(Micro_nelx,Micro_nely,Macro_nelx,Macro_nely,Macro_edofMat_mech,iK_mech,jK_mech,F,U,freedofs_mech, ...
freedofs_Heat,Macro_edofMat_Heat,T,P,iK,jK,Macro_volfrac,Micro_volfrac,Micro_rmin,Macro_rmin,dx,dy,penal,maxIter,option)
%% Preparing Design Variables
clear Macro_x Micro_x
[Macro_x,Micro_x]=get_Initials(Macro_volfrac,Micro_volfrac,Micro_nelx,Micro_nely,Macro_nelx,Macro_nely);
Micro_xPhys = Micro_x;Macro_xPhys=Macro_x;
loop = 0;
while loop < maxIter;
loop = loop+1; 
switch(option)
    case(1);
[Q,dQ]=Homogenization_full(Micro_xPhys,1);
[~,~,~,~,KE_Macro_Homogenized_mech] = elementMatVec(dx/2, dy/2,Q,1);
sK_mech = reshape(KE_Macro_Homogenized_mech(:)*(1e-6+Macro_xPhys(:)'.^penal*(1-1e-6)),64*Macro_nelx*Macro_nely,1);
K_mech = sparse(iK_mech,jK_mech,sK_mech); K_mech = (K_mech+K_mech')/2;
U(freedofs_mech,:) = K_mech(freedofs_mech,freedofs_mech)\F(freedofs_mech,:);
ce_mech = reshape(sum((U(Macro_edofMat_mech)*KE_Macro_Homogenized_mech).*U(Macro_edofMat_mech),2),Macro_nely,Macro_nelx);
c= sum(sum((1e-6+Macro_xPhys.^penal*(1-1e-6)).*ce_mech));
Macro_dc_mech = -penal*(1-1e-6)*Macro_xPhys.^(penal-1).*ce_mech;
Micro_dc_mech = zeros(Micro_nely, Micro_nelx);
for i = 1:Micro_nelx*Micro_nely
dQe_mech = [dQ{1,1}(i) dQ{1,2}(i) dQ{1,3}(i);
            dQ{2,1}(i) dQ{2,2}(i) dQ{2,3}(i);
            dQ{3,1}(i) dQ{3,2}(i) dQ{3,3}(i)];
[~,~,~,~,dKE_mech] = elementMatVec(dx/2, dy/2,dQe_mech,1);
dce_mech = reshape(sum((U(Macro_edofMat_mech)*dKE_mech).*U(Macro_edofMat_mech),2),Macro_nely,Macro_nelx);
Micro_dc_mech(i) = -sum(sum((1e-6+Macro_xPhys.^penal*(1-1e-6)).*dce_mech));
end
[Macro_x, Macro_xPhys, Macro_change] = OC_2D(Macro_x, Macro_dc_mech, Macro_volfrac,Macro_rmin, Macro_nelx,Macro_nely);
[Micro_x, Micro_xPhys, Micro_change] = OC_2D(Micro_x, Micro_dc_mech,  Micro_volfrac,Micro_rmin, Micro_nelx,Micro_nely);

case(2)
[Q,dQ]=Homogenization_full(Micro_xPhys,2);
[~,~,~,~,KE_Homogenized_Heat] = elementMatVec(dx/2, dy/2,Q,2);
sK = reshape(KE_Homogenized_Heat(:)*(1e-6+Macro_xPhys(:)'.^penal*(1-1e-6)),16*Macro_nelx*Macro_nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
T(freedofs_Heat,:) = K(freedofs_Heat,freedofs_Heat)\P(freedofs_Heat,:);
ce = reshape(sum((T(Macro_edofMat_Heat)*KE_Homogenized_Heat).*T(Macro_edofMat_Heat),2),Macro_nely,Macro_nelx);
c= sum(sum((1e-6+Macro_xPhys.^penal*(1-1e-6)).*ce));
Macro_dcH =-penal*(1-1e-6)*Macro_xPhys.^(penal-1).*ce;
Macro_dv = ones(Macro_nely, Macro_nelx);
Micro_dcH = zeros(Micro_nely, Micro_nelx);
for i = 1:Micro_nelx*Micro_nely;  
dQe = [dQ{1,1}(i) dQ{1,2}(i);
dQ{2,1}(i) dQ{2,2}(i)];
[~,~,~,~,dKE] = elementMatVec(dx/2, dy/2,dQe,2);
dce = reshape(sum((T(Macro_edofMat_Heat)*dKE).*T(Macro_edofMat_Heat),2),Macro_nely,Macro_nelx);
Micro_dcH(i) = -sum(sum((1e-6+Macro_xPhys.^penal*(1-1e-6)).*dce));
end
[Macro_x, Macro_xPhys, Macro_change] = OC_2D(Macro_x, Macro_dcH, Macro_volfrac,Macro_rmin, Macro_nelx,Macro_nely);
[Micro_x, Micro_xPhys, Micro_change] = OC_2D(Micro_x, Micro_dcH, Micro_volfrac,Micro_rmin, Micro_nelx,Micro_nely);

end
Macro_xPhys = reshape(Macro_xPhys, Macro_nely, Macro_nelx); Micro_xPhys = reshape(Micro_xPhys, Micro_nely, Micro_nelx);
%% To be deleted when submitted the Program
clf;
%% PRINT RESULTS
fprintf(' It.:%5i Obj.:%11.4f Macro_Vol.:%7.3f Micro_Vol.:%7.3f Macro_ch.:%7.3f Micro_ch.:%7.3f\n',loop,mean(Macro_xPhys(:))...
,mean(Micro_xPhys(:)), Macro_change, Micro_change);
hold on;
switch(option)
case(1);
figure (1);clf;colormap(gray); imagesc(1-Micro_x); caxis([0 1]); 
axis equal; axis off; t1=title('Microscale Design Mechanical');t1.Color ='#E50112 ';drawnow;
figure (2);clf;colormap(gray); imagesc(1-Macro_x); caxis([0 1]);
axis equal; axis off; t2=title('Macroscale Design Mechanical');t2.Color ='#0C95D1';drawnow;
case(2);
figure (4);clf;colormap(gray); imagesc(1-Micro_x); caxis([0 1]); 
axis equal; axis off; t1=title('Microscale Design Heat');t1.Color ='#E600E2';drawnow;
figure (5);clf;colormap(gray); imagesc(1-Macro_x); caxis([0 1]);
axis equal; axis off; t2=title('Macroscale Design Heat');t2.Color ='#09D61F';drawnow;
end;
end;
C_concurrent=c;CH=Q;
disp('********End Concurrent Solver case *********\n');
end