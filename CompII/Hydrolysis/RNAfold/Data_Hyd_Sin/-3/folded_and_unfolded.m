non_folded=0;
folded=0;
for i=1:tf
    disp(i);
    Lt_cell=Lt{i};
    if size(Lt_cell,2)~=0
        for a=1:size(Lt,2)
            cmd  = ['echo ', Lt_cell{a}, ' | RNAfold'];
            [status, result] = system(cmd);
            clean_result = regexprep(result, '\x1B\[[0-9;]*[A-Za-z]', '');
            lines = strsplit(strtrim(clean_result), '\n');
            struct_energy=split(lines{2});
            struct=struct_energy{1};
            energy=struct_energy{3};
            energy=str2num(energy(1:end-1));
            if any(struct~='.')
                folded = folded + 1;
            else
                non_folded = non_folded + 1;
            end
        end
    end
end