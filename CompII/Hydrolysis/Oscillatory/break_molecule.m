function [new_molecules]=break_molecule(breaking_positions, m_break, m_break_OCR, m_break_bl, m_break_E)
    
    new_molecules={};       
    is_mbreak_dsRNA = 0;    
    
    if any(m_break=='_') 
        is_mbreak_dsRNA = 1;
        mark=find(m_break=='_');
        m_left=m_break(1:(mark-1)); % Left strand 
        ml_OCR=m_break_OCR(1:(mark-1)); 
        ml_bl=m_break_bl(1:(mark-1)); 
        ml_E=m_break_E(1:(mark-1)); 
        size_m_break=size(m_left,2); % Size of the strand that is going to break (just to store this information later)

        m_right=m_break((mark+1):end); % Right strand
        mr_OCR=m_break_OCR((mark+1):end);
        mr_E=m_break_E((mark+1):end);
        
        new_bp=[0 breaking_positions, size(m_left,2)];
        
        for h_site=1:(size(new_bp,2)-1)
            
            new_ml=m_left(new_bp(h_site)+1:new_bp(h_site+1));
            new_ml_OCR=ml_OCR(new_bp(h_site)+1:new_bp(h_site+1)); 
            new_ml_bl=ml_bl(new_bp(h_site)+1:new_bp(h_site+1)); % Remember bl labels are stored only for left strand in dsRNA.
            new_ml_bl(end)='X';
            new_ml_E=ml_E(new_bp(h_site)+1:new_bp(h_site+1));
            size_new_ml = size(new_ml,2); % Size of the new molecule that is going to form as the size of its left strand if it is dsRNA
            
            new_mr=m_right(new_bp(h_site)+1:new_bp(h_site+1));
            new_mr_OCR=mr_OCR(new_bp(h_site)+1:new_bp(h_site+1));
            new_mr_E=mr_E(new_bp(h_site)+1:new_bp(h_site+1));
            
            % Construct the new molecules after the break (hydrolysis) m1 and m2:

            if any(new_mr~='-') % If the new molecule m1 is still a dsRNA (for example, it can happen if this molecule 'aaa_-uu' is hydrolyzed in the first p-bond of the left molecule).
                new_m = strcat(new_ml,'_',new_mr);
                new_m_OCR = strcat(new_ml_OCR,'_',new_mr_OCR);
                new_m_bl = new_ml_bl; 
                new_m_E = strcat(new_ml_E,'_',new_mr_E);
                is_m1_dsRNA = 1;
            else % If the new molecule m1 is a ssRNA
                new_m = new_ml;
                new_m_OCR = new_ml_OCR;
                new_m_bl = new_ml_bl;
                new_m_E = new_ml_E;
            end
            
            if size(new_m,2)==1
                new_m_E='!';
                new_m_OCR='O';
            end
            
            new_molecules=[new_molecules {new_m; new_m_OCR; new_m_bl; new_m_E}];
        end
    else
        new_bp=[0 breaking_positions size(m_break,2)];
        for h_site=1:(size(new_bp,2)-1)
            new_m = m_break(new_bp(h_site)+1:new_bp(h_site+1));
            new_m_OCR = m_break_OCR(new_bp(h_site)+1:new_bp(h_site+1));
            new_m_bl = m_break_bl(new_bp(h_site)+1:new_bp(h_site+1));
            new_m_bl(end)='X';
            new_m_E = m_break_E(new_bp(h_site)+1:new_bp(h_site+1));
            size_new_m=size(new_m,2); % Size of the new molecule that forms after hydrolysis
            
            if size(new_m,2)==1
                new_m_E='!';
                new_m_OCR='O';
            end
            
            new_molecules=[new_molecules {new_m; new_m_OCR; new_m_bl; new_m_E}];
        end
    end
    
