function [new_R1_E]=E_table_E1_V2(R1_E,R2_E)

R1_R2=[R1_E; R2_E]; % Remember row 2 is the template.
new_R1_E='';

for i=1:size(R1_R2,2)
    
    if R1_R2(2,i)=='X' % This if is probably unneeded because I believe none of E labels is 'X'
        new_R1_E(i)='X';

    else
        new_R1_E(i)=R1_R2(2,i);
    end
    
end
