function [Macro_x,Micro_x]=get_Initials(Macro_volfrac,Micro_volfrac,Micro_nelx,...
Micro_nely,Macro_nelx,Macro_nely)
Macro_x = repmat(Macro_volfrac,Macro_nely,Macro_nelx);
Micro_x = repmat(Micro_volfrac,Micro_nely,Micro_nelx);
for i = 1:Micro_nelx;for j = 1:Micro_nely;
if sqrt((i-Micro_nelx/2)^2+(j-Micro_nely/2)^2)< min(Micro_nelx,Micro_nely)/4;
Micro_x(j,i) = Micro_volfrac/4;
end;end;end;end