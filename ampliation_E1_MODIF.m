function [Lnew]=ampliation_E1_MODIF(L,R1,OCR_R1,bl_R1,E_R1)

s=size(L,2);
Lnew=cell(3,s+1);
Lnew{1,1}=R1;
Lnew{2,1}=OCR_R1;
Lnew{3,1}=bl_R1;
Lnew{4,1}=E_R1;

Lnew(:,2:(s+1))=L;