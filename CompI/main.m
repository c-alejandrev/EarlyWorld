
function [Lt,Lt_OCR,Lt_bl,tf,Ut,Ut_OCR,Ct,Ct_OCR,Ct_bl,Ot_C,Ot_U,Activity,R_CELL,C_CELL]=TemporalBeta_E1(A,D,N,positions,t,beta,alpha) 
% -- Here beta and alpha are vectors, each position been a different value for a particular timestep

%% INITIALIZING THE VARIABLES

Lambda=positions; %Returns the length of the templante chain C, that is always going to be the total number of positions (e.g. 200)

L=cell(3,sum(N)); % Create a cell, that is a sort of matrix like structure than can harbor any data type. Here they are all strings. 
                  % Here we are creating this cell with just 3 row sand
                  % sum(N) columns, where sum(N) is the total amount of nts
                  % in the media. 
                  % Row 1 => The nt sequence of the ssRNAs or dsRNAs (dsRNAs' chains are separated by '_'). E.g {a} {u} {au_u} {cgcg}
                  % Row 2 => The OCR label of each nt (also dsRNAs' chains are separated by '_'). E.g {O} {O} {OO_C} {CCCC} % Notice that there are O (letter), not 0 (number)
                  % Row 3 => The "bad linkages" (also dsRNAs' chains are separated by '_'). E.g {X} {X} {0X_X} {010X} % The X means that there isn't a link, also, the sequence
                  % starts at the first link, but it has the same length as the sequence of nts, therefore, last element will always be an X.  
                  % Also, here, we use numbers 0 and 1 and the letter X: 0 is a good link and 1 is a bad link.


Lt=''; % Lt is the variable that contains the information of what is the composition of the solution (in terms of monomers, dimers, etc.), i.e., sequences of nts. It is row 1 of L variable.
Lt_OCR=''; % String variable that contains the information of OCR labels. They will match exactly with the sequences of nts in variable Lt. It is row 2 of L variable.
Lt_bl=''; % String variable that contains the information of bad links. They will match exactly with the sequences of nts in variable Lt. It is row 3 of L variable.

Ot_C=[];  % Occupancy of level 0, the clay (C)
Ot_U=[];  % Occupancy of level 1, templates or united nucleosides (U)


Ct=''; % Ct stores the sequence of clay-template linkages at every time step, level 0 across time => E.g. '----...augcc...gccu...---' %Where - means that such position is empty.
Ct_OCR=''; % String that stores the OCR label of the sequences attached to clay at each time step => E.g. '----...OOOOO...CCCC...---'
Ct_bl=''; % String that stores the bad links of the sequences attached to clay at each time step => E.g.  '----...0000X...010X...---' % 1 Means a bad link and X, the absence of link, because there is always one less link than the number of nts
Ut=''; % Ut stores the sequence of template-complementary linkages at every time step, level 1 across time => E.g. '----...uacgg...cg--...---'
Ut_OCR=''; % String that stores the OCR label of the sequences attached to clay at each time step => E.g. '----...CCCCC...RR--...---'

% -- Next loop is just to fill the L cell with all the nucleotides at that are in the solution at the begining of the simulation
    % L structure is more ore less like: {A} {A} ... {C} {C} ... {G} {G} ... {U} {U}
l=1;
for i=1:size(N,2) % For each type of nucleotide (4 normally)
    for n=1:N(i) % For each of the nucleotides of the same kind (2*template lenght normally)
        L{1,l}=A(i);
        L{2,l}='O';
        L{3,l}='X';
        l=l+1;
    end
end

C=''; C_OCR='';C_bl=''; % Remember Lambda is the length of the total number of positions in clay C
C(1:Lambda)='-'; C_OCR(1:Lambda)='-'; C_bl(1:Lambda)='-';

U=''; U_OCR='';
U(1:Lambda)='-';U_OCR(1:Lambda)='-';

O_C=zeros(1,Lambda); % Occupancy of clay (level 0, C) in the form of a matrix of size (1 row, Lambda columns)
O_U=zeros(1,Lambda); % Occupancy of level 1 (U) in the form of a matrix of size (1, Lambda)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables created to measure the activity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Activity=cell(5,3); %Cell to store activity data in the following way
%{                    
                                            All structures        Split by polymer size                               
                      
                               'fixed level 0'             150                [140, 7, 3, 0,...]    % Array is ordered by polymer size 1, 2, 3,...,L      
                              'stripped level 0'           700                [690, 10, 0, 0,...] 
                                'fixed level 1'             150                [140, 7, 3, 0,...]    
                              'stripped level 1'           700                [690, 10, 0, 0,...] 
                              'stripped level 2'           700                [690, 10, 0, 0,...] 
                      %}

total_fixed_0=0;
total_fixed_1=0;
total_stripped_0=0;
total_stripped_1=0;
total_stripped_2=0;
array_fixed_0=zeros(1,Lambda);
array_fixed_1=zeros(1,Lambda);
array_stripped_0=zeros(1,Lambda);
array_stripped_1=zeros(1,Lambda);
array_stripped_2=zeros(1,Lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables created to measure the REPLICANTS AND COMPLEMENTARY CHAINS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_CELL={};
C_CELL={};

counter_rand=0;

%% MAIN CODE: Run a simulation over time

for i=1:t
    %disp(i);
    is_dsRNA=false;        
        
        %% OBTAIN INITIAL IMPORTANT VALUES: Lsize, W, P1, P2, R1, Pn and Pi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Lsize=size(L,2); % Number of colums in L (that is the number of structures that are in the pool at that time)
        W=zeros(1,Lsize); % Weigth of chains in the liquid, because dimers have more probability 
                          % to fall in the template than monomers, so there weight is greater, 
                          % and that of trimers is greater than that of dimers,
                          % etc. And that of dsRNA > ssRNA.

        for ele=1:Lsize
            % At time 1 all are monomers, so they should all weight the same: 1, but then weights can change        
            % REmember that in row 1 the sequences are 
            if any(L{1,ele}=='_') % If there is a dsRNA
                subs=size(find(L{1,ele}=='-'),2);
                W(ele)=size(L{1,ele},2)-1-subs; % We substract the '_' if it is a dsRNA and also the possible empty spaces ('-'). E.g. L{1,ele}='gcg_-g-', then W(ele) = 7-1-2 = 4 nts.
            else
                W(ele)=size(L{1,ele},2);
            end

        end

        Wn=W./sum(W); % Normalization
    %     reset(RandStream.getGlobalStream,sum(100*clock)) 
        P1=sum(rand >= cumsum([0,Wn])); % Random selection according to the weigth of the sample
        P2=randi(Lambda); % Random selected position where the nt or dimer or whatever is going to fall


        R1_L=L{1,P1}; % Nucleotide or chain that falls in to the sample (can be ssRNA or dsRNA)

        if any(R1_L=='_') % If we have dsRNA, choose the strand that is going to fall into clay => Always the one that was facing the clay (the molecule cannot turn around).
            sep=find(R1_L=='_');
            R1=R1_L(1:(sep-1)); % Left strand is always going to be the one that was in direct touch with clay
            R1_U=R1_L((sep+1):end); % Extract also the sequence of the attached complementary strand (right one).

            % Extract the OCR label for both strands
            OCR_both=L{2,P1};
            OCR_left=OCR_both(1:(sep-1));
            OCR_right=OCR_both((sep+1):end);
            OCR=OCR_left;

            is_dsRNA=true;        
        else
            R1=R1_L;
            OCR=L{2,P1};
        end

        % We choose randomly which of the nucleosides of R1 interact 
        % (Pn= position of nucleoside of R1 that interacts,
        %  Pi= position of clay or template with which R1 interacts)

        if size(R1,2)==1 % Size 1 nucleoside
        Pi=P2; % place of interaction  
        Pn=1; 
        else % Size more than 1 nucleoside
        Pn=randi(size(R1,2)); %Pn is the position (nt) of the R1 chain that interacts with the position P2 
        Pi=(P2-Pn+1):(P2+size(R1,2)-Pn); % place of interaction, it is a vector that indicates 
                                         % from and until which position (of clay or template) R1 interacts
        end

        neg=any(Pi<=0); % the position is out in the left edge
        % (P2+size(R1,2)-Pn) <= size(O,2) <=== is the condition for the right
        % edge, for the oligomer R1 to not fall out of the the available positions
        
        
        %% EVALUATE STRUCTURE POSSIBLE FIXATION TO LEVEL 0 OR 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (P2+size(R1,2)-Pn) <= size(O_C,2) & neg==0 % conditions in both edges => only if they are in the total available space (Lambda) they can be joined in that time

            if sum(O_C(Pi))==0 % It falls in clay (level 0) with available place

                Rup=rupture(alpha(i),size(R1,2));
                counter_rand=counter_rand+1;

                if Rup==0 % If links survive to one rupture round, actualize the variables

                    O_C(Pi)=ones(size(Pi));
                    C(Pi)=R1;
                    C_OCR(Pi)=OCR;
                    C_bl(Pi)=L{3,P1};

                    if Pi(1)~=1 % Evaluate the link between R1 and the structure to its left, except if R1 exactly falls at position 1
                        if C_bl(Pi(1)-1)~='-'
                            if C_OCR(Pi(1))=='C' & C_OCR(Pi(1)-1)=='C'
                                C_bl(Pi(1)-1)='1';
                            else
                                C_bl(Pi(1)-1)='0';
                            end
                        end
                    end

                    if Pi(end)~=Lambda % Evaluate the link between R1 and the structure to its right, except if R1 exactly falls at last position (Lambda)
                        if C_bl(Pi(end)+1)~='-'
                            if C_OCR(Pi(end))=='C' & C_OCR(Pi(end)+1)=='C'
                                C_bl(Pi(end))='1';
                            else
                                C_bl(Pi(end))='0';
                            end
                        end
                    end

                    if is_dsRNA % If it is a dsRNA we also have to update U and U_ORC and O_U in addition to O_C, C, C_OCR and C_bl

                        U(Pi)=R1_U;
                        Pi_U=Pi(R1_U~='-'); % Occupied place by R1_U, because it can be an oligo that is smaller than its opposite strand. E.g R1_L='gcg_-g-'. Here if Pi is for example [4 5 6], Pi_U would be just [5]
                        O_U(Pi_U)=ones(size(Pi_U));                
                        U_OCR(Pi)=OCR_right;

                    end            

                    % Delete the structure from the liquid    
                    L=reduction_E1(L,P1); 
                    total_fixed_0=total_fixed_0+1;
                    array_fixed_0(size(R1,2))=array_fixed_0(size(R1,2))+1;

                end

            elseif sum(O_C(Pi))==size(R1,2) % If all positions of level 0 are occupied, but none of level 1 is occupied, R1 can fall in level 1 only if it is ssRNA.

                if sum(O_U(Pi))==0 & is_dsRNA==false % Necessary conditions: available place at level 1 and only ssRNA can attach in level 1

                    % 1st) EVALUATE COMPLEMENTARITY

                    R2=C(Pi); % Chain in clay against with we should evaluate the complementarity of R1
                    RC=complementary2(D,R1,R2); 

                    % 2nd) EVALUATE Prup

                    if RC==1

                         Rup = rupture(beta(i),size(R1,2));
                         counter_rand=counter_rand+1;

                         if Rup==0

                             U(Pi)=R1;
                             O_U(Pi)=ones(size(Pi));

                             % 3rd) CHANGE THE OCR LABELS OF R1 

                             R1_OCR=L{2,P1}; % OCR labels of the chain R1 that fell into template
                             R2_OCR=C_OCR(Pi); % OCR labels of the template R2 into which R1 fell
                             new_R1_OCR=OCR_table_E1(R1_OCR,R2_OCR); % New labels for R1 that emerge from our model table, depending on the complementary labels of R2
                             U_OCR(Pi)=new_R1_OCR; 

                             % Delete the structure from the liquid
                             L=reduction_E1(L,P1);
                             total_fixed_1=total_fixed_1+1;
                             array_fixed_1(size(R1,2))=array_fixed_1(size(R1,2))+1;

                         end
                    end



                end

            else 
                % In any other case (some positions in clay are occupied while others are not), we consider that the attached attempt has failed
            end


        end


        %% EVALUATE RUPTURE PROBABILITIES OF ALL THE STRUCTURES PRESENT IN LEVEL 0, 1 OR POOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -- The order doesn't matter, so we are going to start testing level 1 and 2 (pool) linkages (hydrogen bonds)
        
        % LEVEL 2 (POOL) 
        for pos=1:size(L,2)
            
            if any(L{1,pos}=='_') % If we have dsRNA, evaluate Prup 
        
                ds_oligo=L{1,pos};
                ds_OCR=L{2,pos};
                bl_l=L{3,pos}; % Cause we only store the bad links of level 0
                separ=find(ds_oligo=='_');
                oligo_r=ds_oligo((separ+1):end); % level 1 oligo (to the right of '_' separator)
                OCR_r=ds_OCR((separ+1):end);                
                oligo_l=ds_oligo(1:(separ-1)); % level 0 oligo (to the left of '_' separator)
                OCR_l=ds_OCR(1:(separ-1));

                empty_pos=find(oligo_r=='-');
                empty_pos=[0,empty_pos,length(oligo_r)+1];
                diff_empty_pos=diff(empty_pos)-1;
                reduc=false;

                if size(find(diff_empty_pos>0),2)>1 %If we have for example: oligo_r = '---auuagg----cgaguag-'
                                        
                    % --- 1st) Add new entry for the right ssRNA products that are the other products (in plural) of dsRNA rupture
                    
                    temp=find(oligo_r=='-');
                    temps=size(temp,2);
                    k=zeros(1,temps+1);
                    k(1:temps)=temp;                    

                    if temps>1
                        for z=1:temps
                            if z==1 % First position where there is not a link => First oligo ends here.
                                if k(z)>1

                                    Rup = rupture(beta(i),k(1)-1); 
                                    counter_rand=counter_rand+1;

                                    if Rup==1                            
                                        noligo_r=oligo_r(1:k(1)-1); % Delete possible empty ('-') positions of right oligo
                                        nOCR_r=OCR_r(1:k(1)-1);
                                        if any(nOCR_r=='R')
                                            R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1)};                        
                                        end

                                        if any(OCR_r=='C')
                                            C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1)};
                                        end

                                        bl_r='';
                                        if size(noligo_r,2)==1
                                            bl_r='X';
                                        else
                                            bl_r(1:size(noligo_r,2))='0';
                                            bl_r(end)='X';
                                        end
                                        new_OCR_r=nOCR_r;
                                        new_OCR_r(nOCR_r=='R')='O';
                                        L=ampliation_E1(L,noligo_r,new_OCR_r,bl_r);
                                        oligo_r(1:k(1)-1)='-';
                                        OCR_r(1:k(1)-1)='-';
                                        total_stripped_2=total_stripped_2+1;
                                        array_stripped_2(k(1)-1)=array_stripped_2(k(1)-1)+1; 
                                        reduc=true;
                                    end

                                end

                                if (k(z+1)-k(z))>1 % Evaluate second oligo, because it is not evaluated in following loops
                                    Rup = rupture(beta(i),k(z+1)-k(z)-1); 
                                    counter_rand=counter_rand+1;
                                    if Rup==1 %If there was a break, update the molecules in the solution L

                                        noligo_r=oligo_r((k(z)+1):(k(z+1)-1)); % Delete possible empty ('-') positions of right oligo
                                        nOCR_r=OCR_r((k(z)+1):(k(z+1)-1));
                                        if any(nOCR_r=='R')
                                            R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1))};                        
                                        end

                                        if any(OCR_r=='C')
                                            C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1))};
                                        end

                                        bl_r='';
                                        if size(noligo_r,2)==1
                                            bl_r='X';
                                        else
                                            bl_r(1:size(noligo_r,2))='0';
                                            bl_r(end)='X';
                                        end
                                        new_OCR_r=nOCR_r;
                                        new_OCR_r(nOCR_r=='R')='O';
                                        L=ampliation_E1(L,noligo_r,new_OCR_r,bl_r);
                                        oligo_r((k(z)+1):(k(z+1)-1))='-';
                                        OCR_r((k(z)+1):(k(z+1)-1))='-';
                                        total_stripped_2=total_stripped_2+1;
                                        array_stripped_2(k(z+1)-k(z)-1)=array_stripped_2(k(z+1)-k(z)-1)+1; 
                                        reduc=true;
                                    end
                                end

                            elseif z==temps
                                if k(temps)<size(oligo_l,2)
                                    Rup = rupture(beta(i),size(oligo_l,2)-k(temps)); 
                                    counter_rand=counter_rand+1;
                                    if Rup==1

                                        noligo_r=oligo_r(k(temps)+1:size(oligo_l,2)); % Delete possible empty ('-') positions of right oligo
                                        nOCR_r=OCR_r(k(temps)+1:size(oligo_l,2));
                                        if any(nOCR_r=='R')
                                            R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2))};                        
                                        end

                                        if any(OCR_r=='C')
                                            C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2))};
                                        end

                                        bl_r='';
                                        if size(noligo_r,2)==1
                                            bl_r='X';
                                        else
                                            bl_r(1:size(noligo_r,2))='0';
                                            bl_r(end)='X';
                                        end
                                        new_OCR_r=nOCR_r;
                                        new_OCR_r(nOCR_r=='R')='O';
                                        L=ampliation_E1(L,noligo_r,new_OCR_r,bl_r);
                                        oligo_r(k(temps)+1:size(oligo_l,2))='-';
                                        OCR_r(k(temps)+1:size(oligo_l,2))='-';
                                        total_stripped_2=total_stripped_2+1;
                                        array_stripped_2(size(oligo_l,2)-k(temps))=array_stripped_2(size(oligo_l,2)-k(temps))+1;
                                        reduc=true;
                                    end
                                end

                            else
                                if (k(z+1)-k(z))>1
                                    Rup = rupture(beta(i),k(z+1)-k(z)-1); 
                                    counter_rand=counter_rand+1;
                                    if Rup==1 %If there was a break, update the molecules in the solution L
                                        noligo_r=oligo_r((k(z)+1):(k(z+1)-1)); % Delete possible empty ('-') positions of right oligo
                                        nOCR_r=OCR_r((k(z)+1):(k(z+1)-1));
                                        if any(nOCR_r=='R')
                                            R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1))};                        
                                        end

                                        if any(OCR_r=='C')
                                            C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1))};
                                        end

                                        bl_r='';
                                        if size(noligo_r,2)==1
                                            bl_r='X';
                                        else
                                            bl_r(1:size(noligo_r,2))='0';
                                            bl_r(end)='X';
                                        end
                                        new_OCR_r=nOCR_r;
                                        new_OCR_r(nOCR_r=='R')='O';
                                        L=ampliation_E1(L,noligo_r,new_OCR_r,bl_r);
                                        oligo_r((k(z)+1):(k(z+1)-1))='-';
                                        OCR_r((k(z)+1):(k(z+1)-1))='-';
                                        total_stripped_2=total_stripped_2+1;
                                        array_stripped_2(k(z+1)-k(z)-1)=array_stripped_2(k(z+1)-k(z)-1)+1; 
                                        reduc=true;
                                    end
                                end
                            end

                        end

                    elseif temps==1
                        z=1;
                        if k(z)>1

                            Rup = rupture(beta(i),k(1)-1); 
                            counter_rand=counter_rand+1;

                            if Rup==1                            
                                noligo_r=oligo_r(1:k(1)-1); % Delete possible empty ('-') positions of right oligo
                                nOCR_r=OCR_r(1:k(1)-1);
                                if any(nOCR_r=='R')
                                    R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1)};                        
                                end

                                if any(OCR_r=='C')
                                    C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1)};
                                end

                                bl_r='';
                                if size(noligo_r,2)==1
                                    bl_r='X';
                                else
                                    bl_r(1:size(noligo_r,2))='0';
                                    bl_r(end)='X';
                                end
                                new_OCR_r=nOCR_r;
                                new_OCR_r(nOCR_r=='R')='O';
                                L=ampliation_E1(L,noligo_r,new_OCR_r,bl_r);
                                oligo_r(1:k(1)-1)='-';
                                OCR_r(1:k(1)-1)='-';
                                total_stripped_2=total_stripped_2+1;
                                array_stripped_2(k(1)-1)=array_stripped_2(k(1)-1)+1; 
                                reduc=true;
                            end

                        end

                        if k(temps)<size(oligo_l,2) % Evaluate second oligo, because it is not evaluated in following loops
                            Rup = rupture(beta(i),size(oligo_l,2)-k(temps)); 
                            counter_rand=counter_rand+1;
                            if Rup==1

                                noligo_r=oligo_r(k(temps)+1:size(oligo_l,2)); % Delete possible empty ('-') positions of right oligo
                                nOCR_r=OCR_r(k(temps)+1:size(oligo_l,2));
                                if any(nOCR_r=='R')
                                    R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2))};                        
                                end

                                if any(OCR_r=='C')
                                    C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2))};
                                end

                                bl_r='';
                                if size(noligo_r,2)==1
                                    bl_r='X';
                                else
                                    bl_r(1:size(noligo_r,2))='0';
                                    bl_r(end)='X';
                                end
                                new_OCR_r=nOCR_r;
                                new_OCR_r(nOCR_r=='R')='O';
                                L=ampliation_E1(L,noligo_r,new_OCR_r,bl_r);
                                oligo_r(k(temps)+1:size(oligo_l,2))='-';
                                OCR_r(k(temps)+1:size(oligo_l,2))='-';
                                total_stripped_2=total_stripped_2+1;
                                array_stripped_2(size(oligo_l,2)-k(temps))=array_stripped_2(size(oligo_l,2)-k(temps))+1;
                                reduc=true;
                            end
                        end

                    else            
                         %(if temps==0), the entire chain is not necessary to evaluate here, because we already evaluate it
                    end

                    if reduc
                        
                        % --- 2nd) Add new entry for the RNA that is left, that can be just the left ssRNA or another dsRNA...
                        if all(oligo_r=='-')
                            L=ampliation_E1(L,oligo_l,OCR_l,bl_l);
                        else
                            L=ampliation_E1(L,strcat(oligo_l,'_',oligo_r),strcat(OCR_l,'_',OCR_r),bl_l);
                        end

                        % --- 3rd) Delete the entry in L cell of the original dsRNA
                        isdsOligo = cellfun(@(x)isequal(x,ds_oligo),L(1,:));
                        [row,col] = find(isdsOligo); % Find the original dsRNA position, because we have added new entries to L, and positions changed
                        L=reduction_E1(L,col); 

                    end

                else

                    size_rup=size(find(oligo_r~='-'),2); 
                    
                    if size_rup~=0
                        Rup=rupture(beta(i),size_rup);
                        counter_rand=counter_rand+1;
                        
                        if Rup==1 %%% NOTICE: WE HAVE TO COUNT THE 'R' NTS AND TRANSFORM THEM INTO 'O' IF THE BOND BREAKS
                            % ** También tengo que poner entonces los enlaces que representan a ese oligo, que los contaré todos como buenos, porque se han formado a partir de una template. Si hubiese alguno malo, se contaría a la
                            % hora de guardar la información de los diferentes R formados, pero como luego estos se resetean a O no habría problema, serían enlaces "buenos"
                            % ** Por otro lado, sólo evalúo las 'C' labels si provienen del oligo que fue formado por complementariedad, 
                            % porque las otras para poder formar parte del nivel 0 tuvieron antes que pasar por el pool obligatoriamente, donde ya se habrían contado.
                            noligo_r=oligo_r(oligo_r~='-'); % Delete possible empty ('-') positions of right oligo
                            nOCR_r=OCR_r(OCR_r~='-');
                            if any(OCR_r=='R')
                                R_CELL(end+1,:)={OCR_r, oligo_r, bl_l};                        
                            end
                            
                            if any(OCR_r=='C')
                                C_CELL(end+1,:)={OCR_r, oligo_r, bl_l};
                            end
                            
                            % -- Actualize L (the pool)
                            % --- 1st) Delete the entry in L cell of the dsRNA
                            L=reduction_E1(L,pos); 
                            % --- 2nd) Add two new entries for the two ssRNA that are the products dsRNA rupture
                            L=ampliation_E1(L,oligo_l,OCR_l,bl_l);
                            bl_r='';
                            if size(noligo_r,2)==1
                                bl_r='X';
                            else
                                bl_r(1:size(noligo_r,2))='0';
                                bl_r(end)='X';
                            end
                            new_OCR_r=nOCR_r;
                            new_OCR_r(nOCR_r=='R')='O';
                            L=ampliation_E1(L,noligo_r,new_OCR_r,bl_r);
                            
                        end
                    end
                end
                
            end
            
        end
        
        % LEVEL 1 (hydrogen bonds of dsRNA fixed in clay)
        
        temp=find(O_U==0);
        temps=size(temp,2);
        k=zeros(1,temps+1);
        k(1:temps)=temp;

        if temps>1
            for z=1:temps
                if z==1 % First position where there is not a link => First oligo ends here.
                    if k(z)>1
                        
                        Rup = rupture(beta(i),k(1)-1); 
                        counter_rand=counter_rand+1;
                        
                        if Rup==1                            
                            OCR_chain=U_OCR(1:k(1)-1);
                            nt_chain=U(1:k(1)-1);
                            bl_chain=C_bl(1:k(1)-1);
                            bl_rewritten=bl_chain;
                            bl_rewritten(end)='X';
                            if any(OCR_chain=='R')
                                R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten};                        
                            end

                            if any(OCR_chain=='C')
                                C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten};
                            end
                            
                            bl_new='';
                            if size(nt_chain,2)==1
                                bl_new='X';
                            else
                                bl_new(1:size(nt_chain,2))='0'; % All links are "good", because they were generated by template depending polymerization
                                bl_new(end)='X';
                            end
                            new_OCR_chain=OCR_chain;
                            new_OCR_chain(OCR_chain=='R')='O';     
                            
                            %-- Actualize the variables
                            L=ampliation_E1(L,nt_chain,new_OCR_chain,bl_new);
                            U(1:k(1)-1)='-';
                            O_U(1:k(1)-1)=0;
                            U_OCR(1:k(1)-1)='-';
                            total_stripped_1=total_stripped_1+1;
                            array_stripped_1(k(1)-1)=array_stripped_1(k(1)-1)+1; 
                        end
                        
                    end

                    if (k(z+1)-k(z))>1 % Evaluate second oligo, because it is not evaluated in following loops
                        Rup = rupture(beta(i),k(z+1)-k(z)-1);
                        counter_rand=counter_rand+1;
                        if Rup==1 %If there was a break, update the molecules in the solution L

                            OCR_chain=U_OCR((k(z)+1):(k(z+1)-1));
                            nt_chain=U((k(z)+1):(k(z+1)-1));
                            bl_chain=C_bl((k(z)+1):(k(z+1)-1));
                            bl_rewritten=bl_chain;
                            bl_rewritten(end)='X';

                            if any(OCR_chain=='R')
                                R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten};                        
                            end

                            if any(OCR_chain=='C')
                                C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten};
                            end

                            bl_new='';
                            if size(nt_chain,2)==1
                                bl_new='X';
                            else
                                bl_new(1:size(nt_chain,2))='0'; % All links are "good", because they were generated by template depending polymerization
                                bl_new(end)='X';
                            end
                            new_OCR_chain=OCR_chain;
                            new_OCR_chain(OCR_chain=='R')='O';    

                            L=ampliation_E1(L,nt_chain,new_OCR_chain,bl_new);
                            U((k(z)+1):(k(z+1)-1))='-';
                            O_U((k(z)+1):(k(z+1)-1))=0;
                            U_OCR((k(z)+1):(k(z+1)-1))='-';
                            total_stripped_1=total_stripped_1+1;
                            array_stripped_1(k(z+1)-k(z)-1)=array_stripped_1(k(z+1)-k(z)-1)+1;
                        end
                    end

                elseif z==temps
                    if k(temps)<Lambda
                        Rup = rupture(beta(i),Lambda-k(temps)); 
                        counter_rand=counter_rand+1;
                        if Rup==1

                            OCR_chain=U_OCR(k(temps)+1:Lambda);
                            nt_chain=U(k(temps)+1:Lambda);
                            bl_chain=C_bl(k(temps)+1:Lambda);

                            if any(OCR_chain=='R')
                                R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain};                        
                            end

                            if any(OCR_chain=='C')
                                C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain};
                            end

                            bl_new='';
                            if size(nt_chain,2)==1
                                bl_new='X';
                            else
                                bl_new(1:size(nt_chain,2))='0'; % All links are "good", because they were generated by template depending polymerization
                                bl_new(end)='X';
                            end
                            new_OCR_chain=OCR_chain;
                            new_OCR_chain(OCR_chain=='R')='O';    

                            L=ampliation_E1(L,nt_chain,new_OCR_chain,bl_new);
                            U(k(temps)+1:Lambda)='-';
                            O_U(k(temps)+1:Lambda)=0;
                            U_OCR(k(temps)+1:Lambda)='-';
                            total_stripped_1=total_stripped_1+1;
                            array_stripped_1(Lambda-k(temps))=array_stripped_1(Lambda-k(temps))+1;
                        end
                    end

                else
                    if (k(z+1)-k(z))>1
                        Rup = rupture(beta(i),k(z+1)-k(z)-1);
                        counter_rand=counter_rand+1;
                        if Rup==1 %If there was a break, update the molecules in the solution L

                            OCR_chain=U_OCR((k(z)+1):(k(z+1)-1));
                            nt_chain=U((k(z)+1):(k(z+1)-1));
                            bl_chain=C_bl((k(z)+1):(k(z+1)-1));
                            bl_rewritten=bl_chain;
                            bl_rewritten(end)='X';

                            if any(OCR_chain=='R')
                                R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten};                        
                            end

                            if any(OCR_chain=='C')
                                C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten};
                            end

                            bl_new='';
                            if size(nt_chain,2)==1
                                bl_new='X';
                            else
                                bl_new(1:size(nt_chain,2))='0'; % All links are "good", because they were generated by template depending polymerization
                                bl_new(end)='X';
                            end     
                            new_OCR_chain=OCR_chain;
                            new_OCR_chain(OCR_chain=='R')='O';    

                            L=ampliation_E1(L,nt_chain,new_OCR_chain,bl_new);
                            U((k(z)+1):(k(z+1)-1))='-';
                            O_U((k(z)+1):(k(z+1)-1))=0;
                            U_OCR((k(z)+1):(k(z+1)-1))='-';
                            total_stripped_1=total_stripped_1+1;
                            array_stripped_1(k(z+1)-k(z)-1)=array_stripped_1(k(z+1)-k(z)-1)+1;
                        end
                    end
                end

            end

        elseif temps==1
            z=1;

            if k(z)>1
                        
                Rup = rupture(beta(i),k(1)-1); 
                counter_rand=counter_rand+1;
                
                if Rup==1                            
                    OCR_chain=U_OCR(1:k(1)-1);
                    nt_chain=U(1:k(1)-1);
                    bl_chain=C_bl(1:k(1)-1);
                    bl_rewritten=bl_chain;
                    bl_rewritten(end)='X';
                    if any(OCR_chain=='R')
                        R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten};                        
                    end

                    if any(OCR_chain=='C')
                        C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten};
                    end
                    
                    bl_new='';
                    if size(nt_chain,2)==1
                        bl_new='X';
                    else
                        bl_new(1:size(nt_chain,2))='0'; % All links are "good", because they were generated by template depending polymerization
                        bl_new(end)='X';
                    end
                    new_OCR_chain=OCR_chain;
                    new_OCR_chain(OCR_chain=='R')='O';     
                    
                    %-- Actualize the variables
                    L=ampliation_E1(L,nt_chain,new_OCR_chain,bl_new);
                    U(1:k(1)-1)='-';
                    O_U(1:k(1)-1)=0;
                    U_OCR(1:k(1)-1)='-';
                    total_stripped_1=total_stripped_1+1;
                    array_stripped_1(k(1)-1)=array_stripped_1(k(1)-1)+1; 
                end
                
            end

            if k(temps)<Lambda
                Rup = rupture(beta(i),Lambda-k(temps)); 
                counter_rand=counter_rand+1;
                if Rup==1

                    OCR_chain=U_OCR(k(temps)+1:Lambda);
                    nt_chain=U(k(temps)+1:Lambda);
                    bl_chain=C_bl(k(temps)+1:Lambda);

                    if any(OCR_chain=='R')
                        R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain};                        
                    end

                    if any(OCR_chain=='C')
                        C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain};
                    end

                    bl_new='';
                    if size(nt_chain,2)==1
                        bl_new='X';
                    else
                        bl_new(1:size(nt_chain,2))='0'; % All links are "good", because they were generated by template depending polymerization
                        bl_new(end)='X';
                    end
                    new_OCR_chain=OCR_chain;
                    new_OCR_chain(OCR_chain=='R')='O';    

                    L=ampliation_E1(L,nt_chain,new_OCR_chain,bl_new);
                    U(k(temps)+1:Lambda)='-';
                    O_U(k(temps)+1:Lambda)=0;
                    U_OCR(k(temps)+1:Lambda)='-';
                    total_stripped_1=total_stripped_1+1;
                    array_stripped_1(Lambda-k(temps))=array_stripped_1(Lambda-k(temps))+1;
                end
            end

        else %(if temps==0), the entire chain occupying all positions (e.g. 200 positions) can break
            
            Rup = rupture(beta(i),size(O_U,2));
            counter_rand=counter_rand+1;
            if Rup==1 %If there was a break, update the molecules in the solution L

                OCR_chain=U_OCR;
                nt_chain=U;
                bl_chain=C_bl;

                if any(OCR_chain=='R')
                    R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain};                        
                end

                if any(OCR_chain=='C')
                    C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain};
                end

                bl_new='';
                bl_new(1:size(nt_chain,2))='0'; % All links are "good", because they were generated by template depending polymerization
                bl_new(end)='X';
                new_OCR_chain=OCR_chain;
                new_OCR_chain(OCR_chain=='R')='O';   

                L=ampliation_E1(L,nt_chain,new_OCR_chain,bl_new);
                U(1:Lambda)='-';
                O_U(1:Lambda)=0;
                U_OCR(1:Lambda)='-';
                total_stripped_1=total_stripped_1+1;
                array_stripped_1(size(O_U,2))=array_stripped_1(size(O_U,2))+1;
            end
            
        end

        % LEVEL 0 (ionic bonds of RNA-clay)

        temp=find(O_C==0);
        temps=size(temp,2);
        k=zeros(1,temps+1);
        k(1:temps)=temp;

        if temps>1
            for z=1:temps
                if z==1 % First position where there is not a link => First oligo ends here.
                    if k(z)>1 %There is an oligo occupying position 1 until k(z)-1
                        
                        Rup = rupture(alpha(i),k(1)-1); 
                        counter_rand=counter_rand+1;
                        
                        if Rup==1               
                            if any(U_OCR(1:k(1)-1)~='-') % If we have a dsRNA in that position (level 0, but also some positions (or all) in level 1 are occupied)
                                OCR_chain=strcat(C_OCR(1:k(1)-1),'_',U_OCR(1:k(1)-1));
                                nt_chain=strcat(C(1:k(1)-1),'_',U(1:k(1)-1));
                                bl_chain=C_bl(1:k(1)-1);
                                U(1:k(1)-1)='-';
                                O_U(1:k(1)-1)=0;
                                U_OCR(1:k(1)-1)='-';
                            else
                                OCR_chain=C_OCR(1:k(1)-1);
                                nt_chain=C(1:k(1)-1);
                                bl_chain=C_bl(1:k(1)-1);  
                            end
                     
                            
                            %-- Actualize the variables
                            L=ampliation_E1(L,nt_chain,OCR_chain,bl_chain);
                            C(1:k(1)-1)='-';
                            O_C(1:k(1)-1)=0;
                            C_OCR(1:k(1)-1)='-';
                            C_bl(1:k(1)-1)='-';
                            total_stripped_0=total_stripped_0+1;
                            array_stripped_0(k(1)-1)=array_stripped_0(k(1)-1)+1;
                        end
                        
                    end

                    if (k(z+1)-k(z))>1 % Evaluate second oligo, because it is not evaluated in following loops
                        Rup = rupture(alpha(i),k(z+1)-k(z)-1);
                        counter_rand=counter_rand+1;
                        if Rup==1 %If there was a break, update the molecules in the solution L
                            if any(U_OCR((k(z)+1):(k(z+1)-1))~='-') % If we have a dsRNA in that position (level 0, but also some positions (or all) in level 1 are occupied)
                                OCR_chain=strcat(C_OCR((k(z)+1):(k(z+1)-1)),'_',U_OCR((k(z)+1):(k(z+1)-1)));
                                nt_chain=strcat(C((k(z)+1):(k(z+1)-1)),'_',U((k(z)+1):(k(z+1)-1)));
                                bl_chain=C_bl((k(z)+1):(k(z+1)-1));
                                U((k(z)+1):(k(z+1)-1))='-';
                                O_U((k(z)+1):(k(z+1)-1))=0;
                                U_OCR((k(z)+1):(k(z+1)-1))='-';
                            else
                                OCR_chain=C_OCR((k(z)+1):(k(z+1)-1));
                                nt_chain=C((k(z)+1):(k(z+1)-1));
                                bl_chain=C_bl((k(z)+1):(k(z+1)-1));  
                            end
                     
                            
                            %-- Actualize the variables
                            L=ampliation_E1(L,nt_chain,OCR_chain,bl_chain);
                            C((k(z)+1):(k(z+1)-1))='-';
                            O_C((k(z)+1):(k(z+1)-1))=0;
                            C_OCR((k(z)+1):(k(z+1)-1))='-';
                            C_bl((k(z)+1):(k(z+1)-1))='-';
                            total_stripped_0=total_stripped_0+1;
                            array_stripped_0(k(z+1)-k(z)-1)=array_stripped_0(k(z+1)-k(z)-1)+1;
                        end
                    end

                elseif z==temps
                    if k(temps)<Lambda
                        Rup = rupture(alpha(i),Lambda-k(temps)); 
                        counter_rand=counter_rand+1;
                        if Rup==1
                            if any(U_OCR(k(temps)+1:Lambda)~='-') % If we have a dsRNA in that position (level 0, but also some positions (or all) in level 1 are occupied)
                                OCR_chain=strcat(C_OCR(k(temps)+1:Lambda),'_',U_OCR(k(temps)+1:Lambda));
                                nt_chain=strcat(C(k(temps)+1:Lambda),'_',U(k(temps)+1:Lambda));
                                bl_chain=C_bl(k(temps)+1:Lambda);
                                U(k(temps)+1:Lambda)='-';
                                O_U(k(temps)+1:Lambda)=0;
                                U_OCR(k(temps)+1:Lambda)='-';
                            else
                                OCR_chain=C_OCR(k(temps)+1:Lambda);
                                nt_chain=C(k(temps)+1:Lambda);
                                bl_chain=C_bl(k(temps)+1:Lambda);  
                            end
                     
                            
                            %-- Actualize the variables
                            L=ampliation_E1(L,nt_chain,OCR_chain,bl_chain);
                            C(k(temps)+1:Lambda)='-';
                            O_C(k(temps)+1:Lambda)=0;
                            C_OCR(k(temps)+1:Lambda)='-';
                            C_bl(k(temps)+1:Lambda)='-';
                            total_stripped_0=total_stripped_0+1;
                            array_stripped_0(Lambda-k(temps))=array_stripped_0(Lambda-k(temps))+1;
                        end
                    end
                else
                    if (k(z+1)-k(z))>1
                        Rup = rupture(alpha(i),k(z+1)-k(z)-1);
                        counter_rand=counter_rand+1;
                        if Rup==1 %If there was a break, update the molecules in the solution L
                            if any(U_OCR((k(z)+1):(k(z+1)-1))~='-') % If we have a dsRNA in that position (level 0, but also some positions (or all) in level 1 are occupied)
                                OCR_chain=strcat(C_OCR((k(z)+1):(k(z+1)-1)),'_',U_OCR((k(z)+1):(k(z+1)-1)));
                                nt_chain=strcat(C((k(z)+1):(k(z+1)-1)),'_',U((k(z)+1):(k(z+1)-1)));
                                bl_chain=C_bl((k(z)+1):(k(z+1)-1));
                                U((k(z)+1):(k(z+1)-1))='-';
                                O_U((k(z)+1):(k(z+1)-1))=0;
                                U_OCR((k(z)+1):(k(z+1)-1))='-';
                            else
                                OCR_chain=C_OCR((k(z)+1):(k(z+1)-1));
                                nt_chain=C((k(z)+1):(k(z+1)-1));
                                bl_chain=C_bl((k(z)+1):(k(z+1)-1));  
                            end
                     
                            
                            %-- Actualize the variables
                            L=ampliation_E1(L,nt_chain,OCR_chain,bl_chain);
                            C((k(z)+1):(k(z+1)-1))='-';
                            O_C((k(z)+1):(k(z+1)-1))=0;
                            C_OCR((k(z)+1):(k(z+1)-1))='-';
                            C_bl((k(z)+1):(k(z+1)-1))='-';
                            total_stripped_0=total_stripped_0+1;
                            array_stripped_0(k(z+1)-k(z)-1)=array_stripped_0(k(z+1)-k(z)-1)+1;
                        end
                    end
                end

            end

        elseif temps==1
            z=1;
            if k(z)>1 %There is an oligo occupying position 1 until k(z)-1
                        
                Rup = rupture(alpha(i),k(1)-1); 
                counter_rand=counter_rand+1;

                if Rup==1               
                    if any(U_OCR(1:k(1)-1)~='-') % If we have a dsRNA in that position (level 0, but also some positions (or all) in level 1 are occupied)
                        OCR_chain=strcat(C_OCR(1:k(1)-1),'_',U_OCR(1:k(1)-1));
                        nt_chain=strcat(C(1:k(1)-1),'_',U(1:k(1)-1));
                        bl_chain=C_bl(1:k(1)-1);
                        U(1:k(1)-1)='-';
                        O_U(1:k(1)-1)=0;
                        U_OCR(1:k(1)-1)='-';
                    else
                        OCR_chain=C_OCR(1:k(1)-1);
                        nt_chain=C(1:k(1)-1);
                        bl_chain=C_bl(1:k(1)-1);  
                    end
             
                    
                    %-- Actualize the variables
                    L=ampliation_E1(L,nt_chain,OCR_chain,bl_chain);
                    C(1:k(1)-1)='-';
                    O_C(1:k(1)-1)=0;
                    C_OCR(1:k(1)-1)='-';
                    C_bl(1:k(1)-1)='-';
                    total_stripped_0=total_stripped_0+1;
                    array_stripped_0(k(1)-1)=array_stripped_0(k(1)-1)+1;
                end
                
            end

            if k(temps)<Lambda
                Rup = rupture(alpha(i),Lambda-k(temps)); 
                counter_rand=counter_rand+1;
                if Rup==1
                    if any(U_OCR(k(temps)+1:Lambda)~='-') % If we have a dsRNA in that position (level 0, but also some positions (or all) in level 1 are occupied)
                        OCR_chain=strcat(C_OCR(k(temps)+1:Lambda),'_',U_OCR(k(temps)+1:Lambda));
                        nt_chain=strcat(C(k(temps)+1:Lambda),'_',U(k(temps)+1:Lambda));
                        bl_chain=C_bl(k(temps)+1:Lambda);
                        U(k(temps)+1:Lambda)='-';
                        O_U(k(temps)+1:Lambda)=0;
                        U_OCR(k(temps)+1:Lambda)='-';
                    else
                        OCR_chain=C_OCR(k(temps)+1:Lambda);
                        nt_chain=C(k(temps)+1:Lambda);
                        bl_chain=C_bl(k(temps)+1:Lambda);  
                    end
             
                    
                    %-- Actualize the variables
                    L=ampliation_E1(L,nt_chain,OCR_chain,bl_chain);
                    C(k(temps)+1:Lambda)='-';
                    O_C(k(temps)+1:Lambda)=0;
                    C_OCR(k(temps)+1:Lambda)='-';
                    C_bl(k(temps)+1:Lambda)='-';
                    total_stripped_0=total_stripped_0+1;
                    array_stripped_0(Lambda-k(temps))=array_stripped_0(Lambda-k(temps))+1;
                end
            end


        else %if temps==0
            Rup = rupture(alpha(i),size(O_C,2));
            counter_rand=counter_rand+1;
            if Rup==1 %If there was a break, update the molecules in the solution L
                if any(U_OCR~='-') % If we have a dsRNA in that position (level 0, but also some positions (or all) in level 1 are occupied)
                    OCR_chain=strcat(C_OCR,'_',U_OCR);
                    nt_chain=strcat(C,'_',U);
                    bl_chain=C_bl;
                    U(1:Lambda)='-';
                    O_U(1:Lambda)=0;
                    U_OCR(1:Lambda)='-';
                else
                    OCR_chain=C_OCR;
                    nt_chain=C;
                    bl_chain=C_bl;  
                end
         
                
                %-- Actualize the variables
                L=ampliation_E1(L,nt_chain,OCR_chain,bl_chain);
                C(1:Lambda)='-';
                O_C(1:Lambda)=0;
                C_OCR(1:Lambda)='-';
                C_bl(1:Lambda)='-';
                total_stripped_0=total_stripped_0+1;
                array_stripped_0(size(O_C,2))=array_stripped_0(size(O_C,2))+1;
            end
        end

        % Update the pool and the occupancy arrays
        Ltemp='|';
        LOCRtemp='|';
        Lbltemp='|';
        for n=1:size(L,2)
            tempL=[L{1,n},'|'];
            tempLOCR=[L{2,n},'|'];
            tempLbl=[L{3,n},'|'];
            tempL2=[Ltemp,tempL];
            tempLOCR2=[LOCRtemp,tempLOCR];
            tempLbl2=[Lbltemp,tempLbl];
            Ltemp=tempL2; %Ltemp is like this: '|A|A|T|AG|A|AAC|...', where | separates the different oligos/nts.
            LOCRtemp=tempLOCR2; %LOCRtemp is like this: '|O|O|O|CC|O|OCC|...', where | separates the different oligos/nts.
            Lbltemp=tempLbl2; %Lbltemp is like this: '|0|0|0|0X|0|00X|...', where | separates the different oligos/nts.
        end
        delta=2*sum(N)+1-size(Ltemp,2);
        if delta ~= 0  %If delta is 0, it means that we have nothing in the solution L.
            Ltemp(end+1:end+delta)='|';
            LOCRtemp(end+1:end+delta)='|';
            Lbltemp(end+1:end+delta)='|';
        end

        tf=i;
    

% Update the variables I am going to extract from function (RESULTS)
Lt=[Lt;Ltemp];
if i~=1 & size(LOCRtemp,2)~=size(Lt_OCR(end,:),2)
    LOCRtemp(end+1:size(Lt_OCR(end,:),2))='|';
end
if i~=1 & size(Lbltemp,2)~=size(Lt_bl(end,:),2)
    Lbltemp(end+1:size(Lt_bl(end,:),2))='|';
end
Lt_OCR=[Lt_OCR;LOCRtemp];
Lt_bl=[Lt_bl;Lbltemp];
Ot_C=[Ot_C;O_C];
Ot_U=[Ot_U;O_U];
Ct=[Ct;C];
Ct_OCR=[Ct_OCR;C_OCR];
Ct_bl=[Ct_bl;C_bl];
Ut=[Ut;U];
Ut_OCR=[Ut_OCR;U_OCR];


end

%disp(counter_rand);

% Store activity data
Activity{1,1}='fixed level 0';
Activity{2,1}='stripped level 0';
Activity{1,2}=total_fixed_0;
Activity{2,2}=total_stripped_0;
Activity{1,3}=array_fixed_0;
Activity{2,3}=array_stripped_0;
Activity{3,1}='fixed level 1';
Activity{4,1}='stripped level 1';
Activity{3,2}=total_fixed_1;
Activity{4,2}=total_stripped_1;
Activity{3,3}=array_fixed_1;
Activity{4,3}=array_stripped_1;
Activity{5,1}='stripped level 2';
Activity{5,2}=total_stripped_2;
Activity{5,3}=array_stripped_2;