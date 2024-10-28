function [Lnew]=ampliation_E1(L,R1,OCR_R1,bl_R1)
% Add a new molecule (R1) and its related variables (OCR_R1 and bl_R1)
% to the pool matrix (L)
s=size(L,2);
Lnew=cell(3,s+1);
Lnew{1,1}=R1;
Lnew{2,1}=OCR_R1;
Lnew{3,1}=bl_R1;

Lnew(:,2:(s+1))=L;

