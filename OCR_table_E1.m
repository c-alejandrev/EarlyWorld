function [new_R1_OCR]=OCR_table_E1(R1_OCR,R2_OCR)

R1_R2=[R1_OCR; R2_OCR]; % Remember row 2 is the template.
new_R1_OCR='';

for i=1:size(R1_R2,2)
    
    if R1_R2(2,i)=='O'
        if R1_R2(1,i)=='O'
            new_R1_OCR(i)='C';
        elseif R1_R2(1,i)=='C'
            new_R1_OCR(i)='C';
        else
            disp('Error happened, R1 and R2 must only have O or C labels, not others, check R1_R2 matrix below:');
            disp(R1_R2);
            break
        end
        
    elseif R1_R2(2,i)=='C'
        if R1_R2(1,i)=='O'
            new_R1_OCR(i)='R';
        elseif R1_R2(1,i)=='C'
            new_R1_OCR(i)='R';
        else
            disp('Error happened, R1 and R2 must only have O or C labels, not others, check R1_R2 matrix below:');
            disp(R1_R2);
            break
        end
        
    else
        disp('Error happened, R1 and R2 must only have O or C labels, not others, check R1_R2 matrix below:');
        disp(R1_R2);
        break
    end
    
end
