function [Lnew]=ampliation_E1_MODIF(L,R1,OCR_R1,bl_R1,E_R1)
% Add a new molecule (R1) and its related variables (OCR_R1, bl_R1 and E_R1)
% to the pool matrix (L)
s=size(L,2);

%-- Calculate new folding:
if any(R1=='_') | size(R1)==1
    R1_fold={'.',[],[0]};
else
    cmd = ['echo ', R1, ' | RNAfold'];
    [status, result] = system(cmd);
    clean_result = regexprep(result, '\x1B\[[0-9;]*[A-Za-z]', '');
    lines = strsplit(strtrim(clean_result), '\n');
    struct_energy=split(lines{2});
    struct=struct_energy{1};
    energy=struct_energy{3};
    energy=str2num(energy(1:end-1));
    if all(struct=='.')
        R1_fold={struct,[energy],[0]};
    else
        R1_fold={struct,[energy],[1]};
    end
end

Lnew=cell(5,s+1);
Lnew{1,1}=R1;
Lnew{2,1}=OCR_R1;
Lnew{3,1}=bl_R1;
Lnew{4,1}=E_R1;
Lnew{5,1}=R1_fold;

Lnew(:,2:(s+1))=L;
