function [Lt,Lt_OCR,Lt_bl,Lt_E,tf,Ut,Ut_OCR,Ut_E,Ct,Ct_OCR,Ct_bl,Ct_E,Ot_C,Ot_U,Activity,R_CELL,C_CELL,snapshot_found,snapshot_tmax,length_20]=main(A,D,N,positions,t,beta,alpha,polymer_size) 
% -- Here beta and alpha are vectors, each position been a different value for a particular timestep

%% Important considerations before starting the simulation:
%   -- One realization of Comp.II simulates the replication of an initial ssRNA oligomer of length given by input variable "polymer_size"
%   -- In order to keep track of the identity of the oligomers formed in the simulation (i.e., wheter they are products of 
%   template-dependent polymerization or if they are random newly-generated molecules), each nt is going to have an
%   asigned label (called OCR label), with the following meaning:
%       - Label 'O' corresponds to the nts of the original sequence and all the pool nts initialized at the beginning of every realization
%       - Label 'C' corresponds to those nts that have been created by nucleobase complementarity to 'O'
%       - Label 'R' corresponds to those nts that have been created by nucleobase complementarity to 'C'
%       - 'R' nts are rewritten as 'O' when they return to the pool
%       *For more info on these transformation rules see OCR_table_E1.m function
%   -- To identify the copies of the original sequence and to ensure that the sequence of those copies is correct, two new type of varibles are created:
%       (1) Variables []_bl, help us distinguish if two C nts have been joined randomly at the clay level or thanks to template-dependent polymerization (this label is
%       'X' for single nts, '1' if two C nts have been joined randomly and '0' otherwise)
%       (2) Variables []_E, help us know the copies (Cs or Rs) of the original sequence and if they have been reshuffle. At the beginning of the simulation every nt of the original
%       sequence O is given a character that is ordered following the ASCII code (alphabet order). When copies are formed, they inherit the characters of the nts of their
%       template. That way, if the new copy sequence does not follow the alphabetic order, we can know the sequence has changed. The nts of the newly-generated sequences
%       that have nothing to do with the original are given the label '!'.
%   -- Three levels are important for the simulations:
%       - Level 0 or clay => nts, ssRNA and dsRNA can be adsorbed
%       - Level 1 or complementary level => nts can only be placed if they are complementary to the nts in level 0 (dsRNA cannot be placed here)
%       - Level 2 or pool => nts, ssRNA and dsRNA stay in solution

%% INITIALIZING THE VARIABLES

length_20=0;
snapshot_found={};
snapshot_tmax={};

Lambda=positions; % Number of simulated clay positions (usually Lambda=200)

L=cell(4,sum(N)); % Create a cell, that is a matrix like structure than can harbor any data type. Here they are all strings. 
                  % sum(N) is the total amount of nts in the media at the beginning of the simulation. 
                  % Each column contains one molecule and its related labels.
                  % Row 1 => The nt sequence of the ssRNAs or dsRNAs (dsRNA chains are separated by '_'). E.g {a} {u} {auu_ua} {cgcg} <= Example of 4 molecules present in the pool, two are
                  % single nucleotides, the third one is a dsRNA with 3 nt in one strand and 2 nt in the complementary strand, the fourth is a 4-nt long ssRNA.
                  % Row 2 => The OCR label of each nt of the molecule. E.g {O} {O} {OOO_CC} {CCCC} % Notice that O here is a letter, not the number
                  % Row 3 => The bl label of the molecule. E.g {X} {X} {00X_0X} {010X} % The X means that there isn't a link (because the number of links in a molecule is their length - 1),
                  % also, the sequence starts at the first link, but it has the same length as the sequence of nts, therefore, last element will always be an X.
                  % Row 4 => The E lablel of each nt of the molecule. E.g {!} {!} {AB!_AB} {!!JK}


Lt=''; % Created to store the info of row 1 of L variable for each time step.
Lt_OCR=''; % Created to store the info of row 2 of L variable for each time step.
Lt_bl=''; % Created to store the info of row 3 of L variable for each time step.
Lt_E=''; % Created to store the info of row 4 of L variable for each time step.

Ot_C=[];  % Occupancy of the clay: 0 if the position is empty and 1 otherwise. It is a matrix that stores the occupancy along the entire simulation time (each row is a time step).
Ot_U=[];  % Occupancy of the complementary level: 0 if the position is empty and 1 otherwise. It is a matrix that stores the occupancy along the entire simulation time (each row is a time step).


Ct=''; % State of each clay position (occupied by the nts of the alphabet or empty) at every time step => E.g. '----augcc--gccu---', where '-' means that such position is empty.
Ct_OCR=''; % OCR label of each clay position occupied by a nt at every time step => E.g. '----...OOOOO--CCCC...---'
Ct_bl=''; % bl label of each clay position occupied by a nt at every time step => E.g. '----...0000X--010X...---'
Ct_E=''; % E label of each clay position occupied by a nt at every time step => E.g. '----...!!!!!--AB!!...---'
Ut=''; % State of each complementary (level 1) position (occupied by the nts of the alphabet or empty) at every time step => E.g. '----...uacgg--cg--...---'
Ut_OCR=''; % OCR label of each complementary position occupied by a nt at every time step => E.g. '----...CCCCC...RR--...---'
Ut_E=''; % E label of each complementary position occupied by a nt at every time step => E.g. '----...!!!!!--AB!!...---'



% -- Initializing the variables with their values at the beginning of the simulation
l=1;
for y=1:size(N,2) % For each type of nucleotide (4 normally)
    for n=1:N(y) % For each of the nucleotides of the same kind (2*template lenght normally)
        L{1,l}=A(y);
        L{2,l}='O';
        L{3,l}='X';
        L{4,l}='!'; 
        l=l+1;
    end
end

C=''; C_OCR='';C_bl=''; C_E='';
C(1:Lambda)='-'; C_OCR(1:Lambda)='-'; C_bl(1:Lambda)='-'; C_E(1:Lambda)='-';

U=''; U_OCR=''; U_E='';
U(1:Lambda)='-';U_OCR(1:Lambda)='-'; U_E(1:Lambda)='-'; 

O_C=zeros(1,Lambda); 
O_U=zeros(1,Lambda); 

% -- Create random initial polymer and place it in random clay position

% place_first_polymer=randi(Lambda-polymer_size+1);
% ascii='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
% cont=0;
% for nt=place_first_polymer:(place_first_polymer+polymer_size-1) % This sould be place_first_polymer+polymer_size-1, cause if not, the size is S0+1.
%     cont=cont+1;
%     C(nt)=A(randi(size(A,2))); % Create random polymer with alphabet A and asign it all its labels
%     C_OCR(nt)='O';
%     O_C(nt)=1;
%     C_E(nt)=ascii(cont); 
%     if nt==place_first_polymer+polymer_size
%         C_bl(nt)='X';
%     else
%         C_bl(nt)='0';
%     end
% end


% -- Variables created to measure the activity (number of molecules created along time and split by their lengths)

Activity=cell(5,3); % Activity stored in the following way (example)
%{                    
                                                        All structures        Split by polymer size                               
                      
                          'adsorbed to level 0 (clay)'        150                [140, 7, 3, 0,...]    % Array is ordered by polymer size 1, 2, 3,...,L      
                            'detached from level 0'           700                [690, 10, 0, 0,...] 
                            'hybridized to level 1'           150                [140, 7, 3, 0,...]    
                            'denatured from level 1'          700                [690, 10, 0, 0,...]   
                             'denatured at level 2'           700                [690, 10, 0, 0,...]   
                      %}

total_fixed_0=0; % Total number of molecules (single nts and polymers) adsorbed to the clay along the entire simulation
total_fixed_1=0; % Total number of molecules (single nts and polymers) hybridized to level 1 along the entire simulation
total_stripped_0=0; % Total number of molecules detached from the clay along the entire simulation
total_stripped_1=0; % Total number of molecules denatured from level 1 along the entire simulation
total_stripped_2=0; % Total number of dsRNA molecules denatured at level 2 along the entire simulation
array_fixed_0=zeros(1,Lambda); % Number of molecules adsorbed to the clay along the entire simulation split by polymer length
array_fixed_1=zeros(1,Lambda); % Number of molecules adsorbed to level 1 along the entire simulation split by polymer length
array_stripped_0=zeros(1,Lambda); % Number of molecules detached from the clay along the entire simulation split by polymer length
array_stripped_1=zeros(1,Lambda); % Number of molecules denatured from level 1 along the entire simulation split by polymer length
array_stripped_2=zeros(1,Lambda); % Number of dsRNA molecules denatured at level 2 along the entire simulation split by polymer length


% -- IMPORTANT: Variables created to store all complementary (C) and replicate (R) molecules
R_CELL={}; % Every row stores a molecule with at least one nt being of type 'R' and all of its associated labels (OCR, bl and E)
C_CELL={}; % Every row stores a molecule with at least one nt being of type 'C' and all of its associated labels (OCR, bl and E)

counter_rand=0; % Utilitarian variable

%% MAIN CODE: Run a simulation over time
i=1;
while ~length_20 & i<=t
    
        %disp(i);
        is_dsRNA=false;
          
            %% OBTAIN INITIAL IMPORTANT VALUES: Lsize, W, P1, P2, R1, Pn and Pi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Lsize=size(L,2); % Number of colums in L (i.e., the number of molecules that are in the pool at that time)
            W=zeros(1,Lsize); % Weigth of the molecules of the pool, that is proportional to their length.
    
            for ele=1:Lsize
                % At time step 1 they are all monomers, so they should all weight the same: 1 
                if any(L{1,ele}=='_') % If there is a dsRNA
                    subs=size(find(L{1,ele}=='-'),2);
                    W(ele)=size(L{1,ele},2)-1-subs; % We substract the '_' if it is a dsRNA and also the possible empty spaces ('-'). E.g. L{1,ele}='gcg_-g-', then W(ele) = 7-1-2 = 4 nts.
                else
                    W(ele)=size(L{1,ele},2);
                end
    
            end
    
            Wn=W./sum(W); % Normalization
        %     reset(RandStream.getGlobalStream,sum(100*clock)) 
            P1=sum(rand >= cumsum([0,Wn])); % Random selection of one molecule from the pool as a function of their weigths
            P2=randi(Lambda); % Randomly selected position where the nt or polymer is going to fall
    
    
            R1_L=L{1,P1}; % Extracting the Nt or polymer that is going to fall to the clay
    
            if any(R1_L=='_') % If we have a dsRNA, choose the strand that is going to fall to the clay => Always the one that was facing the clay (the molecule cannot turn around).
                sep=find(R1_L=='_');
                R1=R1_L(1:(sep-1)); % Left strand is always going to be the one that was in direct touch with clay
                R1_U=R1_L((sep+1):end); % Extract also the sequence of the complementary strand (right one).
    
                % Extract the OCR label for both strands
                OCR_both=L{2,P1};
                E_both=L{4,P1};
                OCR_left=OCR_both(1:(sep-1));
                E_left=E_both(1:(sep-1));            
                OCR_right=OCR_both((sep+1):end);
                E_right=E_both((sep+1):end);
                OCR=OCR_left;
                E=E_left;
    
    
                is_dsRNA=true;        
            else
                R1=R1_L;
                OCR=L{2,P1};
                E=L{4,P1};
            end
    
            % We choose randomly which of the nts of R1 molecule interacts with the clay:
            % Pn= position of R1 that interacts, Pi= position(s) of clay with which R1 interacts.
    
            if size(R1,2)==1 % The selected molecule is a single nt
            Pi=P2; % Place of interaction  
            Pn=1; 
            else % The selected molecule is an oligomer
            Pn=randi(size(R1,2)); % Pn is the position (nt) of the R1 chain that interacts with the position P2 
            Pi=(P2-Pn+1):(P2+size(R1,2)-Pn); % Complete place of interaction: it is a vector that indicates all the clay positions with which R1 interacts
            end
    
            neg=any(Pi<=0); % the position is out in the left edge
            % (P2+size(R1,2)-Pn) <= size(O_C,2) % <=== is the condition for the right edge
            
            
            %% EVALUATE STRUCTURE POSSIBLE FIXATION TO LEVEL 0 OR 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            if (P2+size(R1,2)-Pn) <= size(O_C,2) & neg==0 % Conditions in both edges: only if they fall inside the available space (Lambda) they can get adsorbed
    
                if sum(O_C(Pi))==0 % It falls in clay (level 0) with available place
    
                    Rup=rupture(alpha(i),size(R1,2)); % Evaluate desorption probability using the helper function rupture.m
                    counter_rand=counter_rand+1;
    
                    if Rup==0 % If links survive to one rupture round, actualize the variables
    
                        O_C(Pi)=ones(size(Pi));
                        C(Pi)=R1;
                        C_OCR(Pi)=OCR;
                        C_bl(Pi)=L{3,P1};
                        C_E(Pi)=E; 
    
                        if Pi(1)~=1 % Evaluate the polymerization of R1 with the nt to its left (if there is one)
                            if C_bl(Pi(1)-1)~='-'
                                if C_OCR(Pi(1))=='C' & C_OCR(Pi(1)-1)=='C'
                                    C_bl(Pi(1)-1)='1';
                                else
                                    C_bl(Pi(1)-1)='0';
                                end
                            end
                        end
    
                        if Pi(end)~=Lambda % Evaluate the polymerization of R1 with the nt to its right (if there is one)
                            if C_bl(Pi(end)+1)~='-'
                                if C_OCR(Pi(end))=='C' & C_OCR(Pi(end)+1)=='C'
                                    C_bl(Pi(end))='1';
                                else
                                    C_bl(Pi(end))='0';
                                end
                            end
                        end
    
                        if is_dsRNA % If it is a dsRNA we also have to update U, U_ORC and O_U in addition to O_C, C, C_OCR and C_bl
    
                            U(Pi)=R1_U;
                            Pi_U=Pi(R1_U~='-'); % Occupied place by R1_U, because it can be an oligo that is smaller than its opposite strand. E.g R1_L='gcg_-g-'. Here if Pi is for example [4 5 6], Pi_U would be just [5]
                            O_U(Pi_U)=ones(size(Pi_U));                
                            U_OCR(Pi)=OCR_right;
                            U_E(Pi)=E_right;
    
                        end            
    
                        % Delete the falling molecule from the pool   
                        L=reduction_E1(L,P1); 
                        total_fixed_0=total_fixed_0+1;
                        array_fixed_0(size(R1,2))=array_fixed_0(size(R1,2))+1;
    
                    end
    
                elseif sum(O_C(Pi))==size(R1,2) % If all positions at level 0 are occupied, but none at level 1 is occupied, R1 can fall in level 1 only if it is ssRNA.
    
                    if sum(O_U(Pi))==0 & is_dsRNA==false % Necessary conditions: available place at level 1 and only ssRNA can hybridize to level 1
    
                        % 1st) EVALUATE COMPLEMENTARITY
    
                        R2=C(Pi); % RNA molecule in clay against which we should evaluate the complementarity of R1
                        RC=complementary2(D,R1,R2); % Evaluate nucleobase complementarity according to rules set on D and using helper funtion complementary2.m
                        
                        if RC==1 % If there is full complementarity
                             
                             % 2nd) EVALUATE DENATURATION PROBABILITY
    
                             Rup = rupture(beta(i),size(R1,2));
                             counter_rand=counter_rand+1;
    
                             if Rup==0 % If links survive to one rupture round, actualize the variables
    
                                 U(Pi)=R1;
                                 O_U(Pi)=ones(size(Pi));
    
                                 % 3rd) CHANGE THE OCR AND E LABELS OF R1  ACCORDING TO RULES IN OCR_table_E1.m and E_table_E1_V2.m
    
                                 R1_OCR=L{2,P1}; % OCR labels of the molecule R1 that fell to template
                                 R1_E=L{4,P1};
                                 R2_OCR=C_OCR(Pi); % OCR labels of the template R2 above which R1 fell
                                 R2_E=C_E(Pi);
                                 new_R1_OCR=OCR_table_E1(R1_OCR,R2_OCR); 
                                 U_OCR(Pi)=new_R1_OCR; 
                                 new_R1_E=E_table_E1_V2(R1_E,R2_E);
                                 U_E(Pi)=new_R1_E;
    
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

                % -- SEE IF THERE IS A 20-NT (OR GREATER) POLYMER ATTACHED TO CLAY
                padded_O_C=[0,O_C,0];
                diffArray = diff(padded_O_C);
                starts = find(diffArray == 1);
                ends = find(diffArray == -1);
                lengths = ends - starts;

                if any(lengths>=polymer_size)
                    disp('FOUND!')

                    % -- Save the state of the system at this specific time step
                    snapshot_found={};
                    snapshot_found{1}=L;
                    snapshot_found{2}=O_C;
                    snapshot_found{3}=C;
                    snapshot_found{4}=C_OCR;
                    snapshot_found{5}=C_bl;
                    snapshot_found{6}=C_E;
                    snapshot_found{7}=O_U;
                    snapshot_found{8}=U;
                    snapshot_found{9}=U_OCR;
                    snapshot_found{10}=U_E;
                    % Calculate activity cell to store it:
                        Activity_found=cell(5,3);
                        Activity_found{1,1}='fixed level 0';
                        Activity_found{2,1}='stripped level 0';
                        Activity_found{1,2}=total_fixed_0;
                        Activity_found{2,2}=total_stripped_0;
                        Activity_found{1,3}=array_fixed_0;
                        Activity_found{2,3}=array_stripped_0;
                        Activity_found{3,1}='fixed level 1';
                        Activity_found{4,1}='stripped level 1';
                        Activity_found{3,2}=total_fixed_1;
                        Activity_found{4,2}=total_stripped_1;
                        Activity_found{3,3}=array_fixed_1;
                        Activity_found{4,3}=array_stripped_1;
                        Activity_found{5,1}='stripped level 2';
                        Activity_found{5,2}=total_stripped_2;
                        Activity_found{5,3}=array_stripped_2;
                    snapshot_found{11}=Activity_found;
                    snapshot_found{12}=C_CELL;
                    snapshot_found{13}=R_CELL;
                    snapshot_found{14}=i;

                    for x = 1:length(starts) % Find the position in O_C where the >= 20-nt sequence is positioned
                        len = ends(x) - starts(x);
                        if len >=polymer_size
                            new_pol_size=len;
                            start_pos = starts(x);
                            end_pos = ends(x)-1;
                        end
                    end

                    ascii='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
                    cont=0;
                    for nt=start_pos:end_pos 
                        cont=cont+1;
                        %C(nt)=A(randi(size(A,2))); % Create random polymer with alphabet A and asign it all its labels
                        C_OCR(nt)='O';
                        O_C(nt)=1;
                        C_E(nt)=ascii(cont); 
                        if nt==end_pos
                            C_bl(nt)='X';
                        else
                            C_bl(nt)='0';
                        end
                        if O_U(nt)==1
                            U_E(nt)=ascii(cont);
                        end
                    end

                    snapshot_found{15}=C_E;
                    snapshot_found{16}=U_E;
                    snapshot_found{17}=C(start_pos:end_pos); % The new original sequence polymer
                    snapshot_found{18}=new_pol_size; % Its size
                    length_20=1;
                    t_found=i;

                    continue

                end
    
    
            end
    
    
            %% EVALUATE RUPTURE PROBABILITIES OF ALL THE MOLECULES IN LEVEL 0, 1 OR 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % -- The order is not important, so we are going to start testing level 1 and 2 (pool) linkages (hydrogen bonds)
            
            % LEVEL 2 (POOL) => Evaluate the denaturation probability of all dsRNA molecules present in the pool in the current time step
            for pos=1:size(L,2)
                
                if any(L{1,pos}=='_') % If we have dsRNA, evaluate denaturation probability
            
                    ds_oligo=L{1,pos};
                    ds_OCR=L{2,pos};
                    ds_E=L{4,pos};
                    bl_l=L{3,pos}; % Cause we only store the bl label of level 0
                    separ=find(ds_oligo=='_');
                    oligo_r=ds_oligo((separ+1):end); % level 1 oligomer (to the right of '_' separator)
                    OCR_r=ds_OCR((separ+1):end);  
                    E_r=ds_E((separ+1):end); 
                    oligo_l=ds_oligo(1:(separ-1)); % level 0 oligomer (to the left of '_' separator)
                    OCR_l=ds_OCR(1:(separ-1));
                    E_l=ds_E(1:(separ-1));
    
                    empty_pos=find(oligo_r=='-');
                    empty_pos=[0,empty_pos,length(oligo_r)+1];
                    diff_empty_pos=diff(empty_pos)-1;
    %                 diff_empty_pos=diff(empty_pos);
                    reduc=false; 
    
                    if size(find(diff_empty_pos>0),2)>1 % If we have for example: oligo_r = '---auuagg----cgaguag-'
                                            
                        % --- 1st) Add new entry for the right ssRNA products that are the products (in plural) of dsRNA rupture
                        
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
                                            nE_r=E_r(1:k(1)-1);
                                            if any(nOCR_r=='R')
                                                R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r};                        
                                            end
    
                                            if any(OCR_r=='C')
                                                C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r};
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
                                            if k(1)-1==1
                                                nE_r='!';
                                            end
                                            L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                            oligo_r(1:k(1)-1)='-';
                                            OCR_r(1:k(1)-1)='-';
                                            E_r(1:k(1)-1)='-';
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
                                            nE_r=E_r((k(z)+1):(k(z+1)-1));
                                            
                                            if any(nOCR_r=='R')
                                                R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r};                        
                                            end
    
                                            if any(OCR_r=='C')
                                                C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r};
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
                                            if (k(z+1)-k(z)-1)==1
                                                nE_r='!';
                                            end
                                            L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                            oligo_r((k(z)+1):(k(z+1)-1))='-';
                                            OCR_r((k(z)+1):(k(z+1)-1))='-';
                                            E_r((k(z)+1):(k(z+1)-1))='-';
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
                                            nE_r=E_r(k(temps)+1:size(oligo_l,2));
    
                                            if any(nOCR_r=='R')
                                                R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r};                        
                                            end
    
                                            if any(OCR_r=='C')
                                                C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r};
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
                                            if (size(oligo_l,2)-k(temps))==1
                                                nE_r='!';
                                            end
                                            L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                            oligo_r(k(temps)+1:size(oligo_l,2))='-';
                                            OCR_r(k(temps)+1:size(oligo_l,2))='-';
                                            E_r(k(temps)+1:size(oligo_l,2))='-';
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
                                            nE_r=E_r((k(z)+1):(k(z+1)-1));
    
                                            if any(nOCR_r=='R')
                                                R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r};                        
                                            end
    
                                            if any(OCR_r=='C')
                                                C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r};
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
                                            if (k(z+1)-k(z)-1)==1
                                                nE_r='!';
                                            end
                                            L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                            oligo_r((k(z)+1):(k(z+1)-1))='-';
                                            OCR_r((k(z)+1):(k(z+1)-1))='-';
                                            E_r((k(z)+1):(k(z+1)-1))='-';
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
                                    nE_r=E_r(1:k(1)-1);
    
                                    if any(nOCR_r=='R')
                                        R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r};                        
                                    end
    
                                    if any(OCR_r=='C')
                                        C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r};
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
                                    if (k(1)-1)==1
                                        nE_r='!';
                                    end
                                    L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                    oligo_r(1:k(1)-1)='-';
                                    OCR_r(1:k(1)-1)='-';
                                    E_r(1:k(1)-1)='-';
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
                                    nE_r=E_r(k(temps)+1:size(oligo_l,2));
    
                                    if any(nOCR_r=='R')
                                        R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r};                        
                                    end
    
                                    if any(OCR_r=='C')
                                        C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r};
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
                                    if (size(oligo_l,2)-k(temps))==1
                                        nE_r='!';
                                    end
                                    L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                    oligo_r(k(temps)+1:size(oligo_l,2))='-';
                                    OCR_r(k(temps)+1:size(oligo_l,2))='-';
                                    E_r(k(temps)+1:size(oligo_l,2))='-';
                                    total_stripped_2=total_stripped_2+1;
                                    array_stripped_2(size(oligo_l,2)-k(temps))=array_stripped_2(size(oligo_l,2)-k(temps))+1;
                                    reduc=true;
                                end
                            end
    
    
                        else
                             %(if temps==0), it is not necessary to evaluate the entire chain rupture here, because we already evaluated it
                        end
    
                        if reduc
                            
                            % --- 2nd) Add new entry for the RNA that is left, that can be just the left ssRNA or another dsRNA
                            if all(oligo_r=='-')
                                if size(E_l,2)==1
                                    E_l='!';
                                end
                                L=ampliation_E1_MODIF(L,oligo_l,OCR_l,bl_l,E_l);
                            else
                                L=ampliation_E1_MODIF(L,strcat(oligo_l,'_',oligo_r),strcat(OCR_l,'_',OCR_r),bl_l,strcat(E_l,'_',E_r));
                            end
    
                            % --- 3rd) Delete the entry in L cell of the original dsRNA
                            isdsOligo = cellfun(@(x)isequal(x,ds_oligo),L(1,:));
                            [row,col] = find(isdsOligo); % Find the original dsRNA position, because we have added new entries to L, and positions changed
                            L=reduction_E1(L,col); 
    
                        end
    
                    else % If oligo_r is simpler
    
                        size_rup=size(find(oligo_r~='-'),2); 
                        
                        if size_rup~=0
                            Rup=rupture(beta(i),size_rup);
                            counter_rand=counter_rand+1;
                            
                            if Rup==1 
                                noligo_r=oligo_r(oligo_r~='-'); % Delete possible empty ('-') positions of right oligo
                                nOCR_r=OCR_r(OCR_r~='-');
                                nE_r=E_r(E_r~='-');
                                if size(nE_r,2)==1
                                    nE_r='!';
                                end
    
                                if any(OCR_r=='R')
                                    R_CELL(end+1,:)={OCR_r, oligo_r, bl_l, E_r};                        
                                end
                                
                                if any(OCR_r=='C')
                                    C_CELL(end+1,:)={OCR_r, oligo_r, bl_l, E_r};
                                end
                                
                                % -- Actualize L (the pool)
                                % --- 1st) Delete the entry in L cell of the dsRNA
                                L=reduction_E1(L,pos); 
                                % --- 2nd) Add two new entries for the two ssRNA that are the products dsRNA rupture
                                if size(E_l,2)==1
                                    E_l='!';
                                end
                                L=ampliation_E1_MODIF(L,oligo_l,OCR_l,bl_l,E_l);
                                bl_r='';
                                if size(noligo_r,2)==1
                                    bl_r='X';
                                else
                                    bl_r(1:size(noligo_r,2))='0';
                                    bl_r(end)='X';
                                end
                                new_OCR_r=nOCR_r;
                                new_OCR_r(nOCR_r=='R')='O';
                                L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                
                            end
                        end
                    end
                    
                end
                
            end
            
            % LEVEL 1 => Evaluate the denaturation probability of the molecules hybridized to level 1 in the current time step
            
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
                                E_chain=U_E(1:k(1)-1);
    
                                if any(OCR_chain=='R')
                                    R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};                        
                                end
    
                                if any(OCR_chain=='C')
                                    C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};
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
                                if (k(1)-1)==1
                                    E_chain='!';
                                end
                                L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                                U(1:k(1)-1)='-';
                                O_U(1:k(1)-1)=0;
                                U_OCR(1:k(1)-1)='-';
                                U_E(1:k(1)-1)='-';
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
                                E_chain=U_E((k(z)+1):(k(z+1)-1));
    
                                if any(OCR_chain=='R')
                                    R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten,E_chain};                        
                                end
    
                                if any(OCR_chain=='C')
                                    C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten,E_chain};
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
    
                                if (k(z+1)-k(z)-1)==1
                                    E_chain='!';
                                end
                                L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                                U((k(z)+1):(k(z+1)-1))='-';
                                O_U((k(z)+1):(k(z+1)-1))=0;
                                U_OCR((k(z)+1):(k(z+1)-1))='-';
                                U_E((k(z)+1):(k(z+1)-1))='-';
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
                                E_chain=U_E(k(temps)+1:Lambda);
    
                                if any(OCR_chain=='R')
                                    R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};                        
                                end
    
                                if any(OCR_chain=='C')
                                    C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};
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
    
                                if (Lambda-k(temps))==1
                                    E_chain='!';
                                end
                                L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                                U(k(temps)+1:Lambda)='-';
                                O_U(k(temps)+1:Lambda)=0;
                                U_OCR(k(temps)+1:Lambda)='-';
                                U_E(k(temps)+1:Lambda)='-';
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
                                E_chain=U_E((k(z)+1):(k(z+1)-1));                           
    
                                if any(OCR_chain=='R')
                                    R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};                        
                                end
    
                                if any(OCR_chain=='C')
                                    C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};
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
    
                                if (k(z+1)-k(z)-1)==1
                                    E_chain='!';
                                end
                                L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                                U((k(z)+1):(k(z+1)-1))='-';
                                O_U((k(z)+1):(k(z+1)-1))=0;
                                U_OCR((k(z)+1):(k(z+1)-1))='-';
                                U_E((k(z)+1):(k(z+1)-1))='-';
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
                        E_chain=U_E(1:k(1)-1);
    
                        if any(OCR_chain=='R')
                            R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};                        
                        end
    
                        if any(OCR_chain=='C')
                            C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};
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
                        if (k(1)-1)==1
                            E_chain='!';
                        end
                        L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                        U(1:k(1)-1)='-';
                        O_U(1:k(1)-1)=0;
                        U_OCR(1:k(1)-1)='-';
                        U_E(1:k(1)-1)='-';
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
                        E_chain=U_E(k(temps)+1:Lambda);
                        if any(OCR_chain=='R')
                            R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};                        
                        end
    
                        if any(OCR_chain=='C')
                            C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};
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
    
                        if (Lambda-k(temps))==1
                            E_chain='!';
                        end
                        L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                        U(k(temps)+1:Lambda)='-';
                        O_U(k(temps)+1:Lambda)=0;
                        U_OCR(k(temps)+1:Lambda)='-';
                        U_E(k(temps)+1:Lambda)='-';
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
                    E_chain=U_E;
                    if any(OCR_chain=='R')
                        R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};                        
                    end
    
                    if any(OCR_chain=='C')
                        C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};
                    end
    
                    bl_new='';
                    bl_new(1:size(nt_chain,2))='0'; % All links are "good", because they were generated by template depending polymerization
                    bl_new(end)='X';
                    new_OCR_chain=OCR_chain;
                    new_OCR_chain(OCR_chain=='R')='O';   
    
                    L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                    U(1:Lambda)='-';
                    O_U(1:Lambda)=0;
                    U_OCR(1:Lambda)='-';
                    U_E(1:Lambda)='-';
                    total_stripped_1=total_stripped_1+1;
                    array_stripped_1(size(O_U,2))=array_stripped_1(size(O_U,2))+1;
                end
                
            end
    
            % LEVEL 0 => Evaluate the desorption probability of all the molecules adsorbed to level 0
    
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
                                    E_chain=strcat(C_E(1:k(1)-1),'_',U_E(1:k(1)-1));
                                    U(1:k(1)-1)='-';
                                    O_U(1:k(1)-1)=0;
                                    U_OCR(1:k(1)-1)='-';
                                    U_E(1:k(1)-1)='-';
                                else
                                    OCR_chain=C_OCR(1:k(1)-1);
                                    nt_chain=C(1:k(1)-1);
                                    bl_chain=C_bl(1:k(1)-1);  
                                    
                                    if (k(1)-1)==1
                                        E_chain='!';
                                    else
                                        E_chain=C_E(1:k(1)-1);
                                    end
                                end
                         
                                
                                %-- Actualize the variables
                                L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                                C(1:k(1)-1)='-';
                                O_C(1:k(1)-1)=0;
                                C_OCR(1:k(1)-1)='-';
                                C_bl(1:k(1)-1)='-';
                                C_E(1:k(1)-1)='-';
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
                                    E_chain=strcat(C_E((k(z)+1):(k(z+1)-1)),'_',U_E((k(z)+1):(k(z+1)-1)));
                                    U((k(z)+1):(k(z+1)-1))='-';
                                    O_U((k(z)+1):(k(z+1)-1))=0;
                                    U_OCR((k(z)+1):(k(z+1)-1))='-';
                                    U_E((k(z)+1):(k(z+1)-1))='-';
                                else
                                    OCR_chain=C_OCR((k(z)+1):(k(z+1)-1));
                                    nt_chain=C((k(z)+1):(k(z+1)-1));
                                    bl_chain=C_bl((k(z)+1):(k(z+1)-1)); 
                                     
                                    if (k(z+1)-k(z)-1)==1
                                        E_chain='!';
                                    else
                                        E_chain=C_E((k(z)+1):(k(z+1)-1));
                                    end
                                end
                         
                                
                                %-- Actualize the variables
                                L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                                C((k(z)+1):(k(z+1)-1))='-';
                                O_C((k(z)+1):(k(z+1)-1))=0;
                                C_OCR((k(z)+1):(k(z+1)-1))='-';
                                C_bl((k(z)+1):(k(z+1)-1))='-';
                                C_E((k(z)+1):(k(z+1)-1))='-';
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
                                    E_chain=strcat(C_E(k(temps)+1:Lambda),'_',U_E(k(temps)+1:Lambda));
                                    U(k(temps)+1:Lambda)='-';
                                    O_U(k(temps)+1:Lambda)=0;
                                    U_OCR(k(temps)+1:Lambda)='-';
                                    U_E(k(temps)+1:Lambda)='-';
                                else
                                    OCR_chain=C_OCR(k(temps)+1:Lambda);
                                    nt_chain=C(k(temps)+1:Lambda);
                                    bl_chain=C_bl(k(temps)+1:Lambda);  
                                    if (Lambda-k(temps))==1
                                        E_chain='!';
                                    else
                                        E_chain=C_E(k(temps)+1:Lambda);  
                                    end
                                end
                         
                                
                                %-- Actualize the variables
                                L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                                C(k(temps)+1:Lambda)='-';
                                O_C(k(temps)+1:Lambda)=0;
                                C_OCR(k(temps)+1:Lambda)='-';
                                C_bl(k(temps)+1:Lambda)='-';
                                U_E(k(temps)+1:Lambda)='-';
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
                                    E_chain=strcat(C_E((k(z)+1):(k(z+1)-1)),'_',U_E((k(z)+1):(k(z+1)-1)));
                                    U((k(z)+1):(k(z+1)-1))='-';
                                    O_U((k(z)+1):(k(z+1)-1))=0;
                                    U_OCR((k(z)+1):(k(z+1)-1))='-';
                                    U_E((k(z)+1):(k(z+1)-1))='-';
                                else
                                    OCR_chain=C_OCR((k(z)+1):(k(z+1)-1));
                                    nt_chain=C((k(z)+1):(k(z+1)-1));
                                    bl_chain=C_bl((k(z)+1):(k(z+1)-1)); 
                                    if (k(z+1)-k(z)-1)==1
                                        E_chain='!';
                                    else
                                        E_chain=C_E((k(z)+1):(k(z+1)-1));
                                    end
                                end
                         
                                
                                %-- Actualize the variables
                                L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                                C((k(z)+1):(k(z+1)-1))='-';
                                O_C((k(z)+1):(k(z+1)-1))=0;
                                C_OCR((k(z)+1):(k(z+1)-1))='-';
                                C_bl((k(z)+1):(k(z+1)-1))='-';
                                C_E((k(z)+1):(k(z+1)-1))='-';
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
                            E_chain=strcat(C_E(1:k(1)-1),'_',U_E(1:k(1)-1));
                            U(1:k(1)-1)='-';
                            O_U(1:k(1)-1)=0;
                            U_OCR(1:k(1)-1)='-';
                            U_E(1:k(1)-1)='-';
                        else
                            OCR_chain=C_OCR(1:k(1)-1);
                            nt_chain=C(1:k(1)-1);
                            bl_chain=C_bl(1:k(1)-1); 
                            if (k(1)-1)==1
                                E_chain='!';
                            else
                                E_chain=C_E(1:k(1)-1);
                            end
                        end
                 
                        
                        %-- Actualize the variables
                        L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                        C(1:k(1)-1)='-';
                        O_C(1:k(1)-1)=0;
                        C_OCR(1:k(1)-1)='-';
                        C_bl(1:k(1)-1)='-';
                        C_E(1:k(1)-1)='-';
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
                            E_chain=strcat(C_E(k(temps)+1:Lambda),'_',U_E(k(temps)+1:Lambda));
                            U(k(temps)+1:Lambda)='-';
                            O_U(k(temps)+1:Lambda)=0;
                            U_OCR(k(temps)+1:Lambda)='-';
                            U_E(k(temps)+1:Lambda)='-';
                        else
                            OCR_chain=C_OCR(k(temps)+1:Lambda);
                            nt_chain=C(k(temps)+1:Lambda);
                            bl_chain=C_bl(k(temps)+1:Lambda);  
                            if (Lambda-k(temps))==1
                                E_chain='!';
                            else
                                E_chain=C_E(k(temps)+1:Lambda);
                            end
                        end
                 
                        
                        %-- Actualize the variables
                        L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                        C(k(temps)+1:Lambda)='-';
                        O_C(k(temps)+1:Lambda)=0;
                        C_OCR(k(temps)+1:Lambda)='-';
                        C_bl(k(temps)+1:Lambda)='-';
                        C_E(k(temps)+1:Lambda)='-';
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
                        E_chain=strcat(C_E,'_',U_E);
                        U(1:Lambda)='-';
                        O_U(1:Lambda)=0;
                        U_OCR(1:Lambda)='-';
                        U_E(1:Lambda)='-';
                    else
                        OCR_chain=C_OCR;
                        nt_chain=C;
                        bl_chain=C_bl; 
                        E_chain=C_E;
                    end
             
                    
                    %-- Actualize the variables
                    L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                    C(1:Lambda)='-';
                    O_C(1:Lambda)=0;
                    C_OCR(1:Lambda)='-';
                    C_bl(1:Lambda)='-';
                    C_E(1:Lambda)='-';
                    total_stripped_0=total_stripped_0+1;
                    array_stripped_0(size(O_C,2))=array_stripped_0(size(O_C,2))+1;
                end
            end
    
   
            tf=i; % Utilitarian variable

            if i==t
                snapshot_tmax={};
                % Calculate activity cell to store it:
                        Activity_tmax=cell(5,3);
                        Activity_tmax{1,1}='fixed level 0';
                        Activity_tmax{2,1}='stripped level 0';
                        Activity_tmax{1,2}=total_fixed_0;
                        Activity_tmax{2,2}=total_stripped_0;
                        Activity_tmax{1,3}=array_fixed_0;
                        Activity_tmax{2,3}=array_stripped_0;
                        Activity_tmax{3,1}='fixed level 1';
                        Activity_tmax{4,1}='stripped level 1';
                        Activity_tmax{3,2}=total_fixed_1;
                        Activity_tmax{4,2}=total_stripped_1;
                        Activity_tmax{3,3}=array_fixed_1;
                        Activity_tmax{4,3}=array_stripped_1;
                        Activity_tmax{5,1}='stripped level 2';
                        Activity_tmax{5,2}=total_stripped_2;
                        Activity_tmax{5,3}=array_stripped_2;
               snapshot_tmax{1}=Activity_tmax;
               snapshot_tmax{2}=C_CELL;
               snapshot_tmax{3}=R_CELL;
               snapshot_tmax{4}=C;
               snapshot_tmax{5}=C_OCR;
               snapshot_tmax{6}=C_bl;
               snapshot_tmax{7}=C_E;
            end
    
    
    i=i+1;
        % ---- Update variables. NOTICE: Uncomment following lines to obtain all possible data from simulations in addition to the Activity, R_CELL and C_CELL variables (althoug notice that this will require additional time per simulation)----
    
    %         % Update the pool and the occupancy arrays
    %         Ltemp='|';
    %         LOCRtemp='|';
    %         Lbltemp='|';
    %         LEtemp='|';
    %         for n=1:size(L,2)
    %             tempL=[L{1,n},'|'];
    %             tempLOCR=[L{2,n},'|'];
    %             tempLbl=[L{3,n},'|'];
    %             tempLE=[L{4,n},'|'];
    %             tempL2=[Ltemp,tempL];
    %             tempLOCR2=[LOCRtemp,tempLOCR];
    %             tempLbl2=[Lbltemp,tempLbl];
    %             tempLE2=[LEtemp,tempLE];
    %             Ltemp=tempL2; %Ltemp has a structure like this: '|A|A|T|AG|A|AAC|...', where | separates the different oligos/nts.
    %             LOCRtemp=tempLOCR2; %LOCRtemp has a structure like this: '|O|O|O|CC|O|OCC|...', where | separates the different oligos/nts.
    %             Lbltemp=tempLbl2; %Lbltemp has a structure like this: '|0|0|0|0X|0|00X|...', where | separates the different oligos/nts.
    %             LEtemp=tempLE2; %Lbltemp has a structure like this: '|!|!|!|AB|!|!JK|...', where | separates the different oligos/nts.
    %         end
    % 
%              tf=i; % Utilitarian variable
    % 
    %     if i~=1
    % 
    %         %% -- Lt
    % 
    %         if size(Ltemp,2)==size(Lt(end,:),2)
    %             Lt=[Lt;Ltemp];
    %         else
    %             prev_size=size(Lt(end,:),2);
    %             new_size=size(Ltemp,2);
    % 
    %             if prev_size>new_size
    %                 Ltemp(end+1:size(Lt(end,:),2))='|';
    %             elseif new_size>prev_size
    %                 Lt(1:end,end+1:size(Ltemp(end,:),2))='|';
    %             end
    %             Lt=[Lt;Ltemp];
    %         end
    % 
    %          %% -- LtOCR
    %         
    %         if size(LOCRtemp,2)==size(Lt_OCR(end,:),2)
    %             Lt_OCR=[Lt_OCR;LOCRtemp];
    %         else
    %             prev_size=size(Lt_OCR(end,:),2);
    %             new_size=size(LOCRtemp,2);
    % 
    %             if prev_size>new_size
    %                 LOCRtemp(end+1:size(Lt_OCR(end,:),2))='|';
    %             elseif new_size>prev_size
    %                 Lt_OCR(1:end,end+1:size(LOCRtemp(end,:),2))='|';
    %             end
    %             Lt_OCR=[Lt_OCR;LOCRtemp];
    %         end
    % 
    %         %% -- Ltbl
    % 
    %         if size(Lbltemp,2)==size(Lt_bl(end,:),2)
    %             Lt_bl=[Lt_bl;Lbltemp];
    %         else
    %             prev_size=size(Lt_bl(end,:),2);
    %             new_size=size(Lbltemp,2);
    % 
    %             if prev_size>new_size
    %                 Lbltemp(end+1:size(Lt_bl(end,:),2))='|';
    %             elseif new_size>prev_size
    %                 Lt_bl(1:end,end+1:size(Lbltemp(end,:),2))='|';
    %             end
    %             Lt_bl=[Lt_bl;Lbltemp];
    %         end
    % 
    %         %% -- LtE
    % 
    %         if size(LEtemp,2)==size(Lt_E(end,:),2)
    %             Lt_E=[Lt_E;LEtemp];
    %         else
    %             prev_size=size(Lt_E(end,:),2);
    %             new_size=size(LEtemp,2);
    % 
    %             if prev_size>new_size
    %                 LEtemp(end+1:size(Lt_E(end,:),2))='|';
    %             elseif new_size>prev_size
    %                 Lt_E(1:end,end+1:size(LEtemp(end,:),2))='|';
    %             end
    %             Lt_E=[Lt_E;LEtemp];
    %         end
    % 
    % 
    %     else
    %         Lt=[Lt;Ltemp];
    %         Lt_OCR=[Lt_OCR;LOCRtemp];
    %         Lt_bl=[Lt_bl;Lbltemp];
    %         Lt_E=[Lt_E;LEtemp];
    % 
    %     end
    
%     Ot_C=[Ot_C;O_C];
%     Ot_U=[Ot_U;O_U];
%     Ct=[Ct;C];
%     Ct_OCR=[Ct_OCR;C_OCR];
    % Ct_bl=[Ct_bl;C_bl];
%     Ct_E=[Ct_E;C_E];
%     Ut=[Ut;U];
%     Ut_OCR=[Ut_OCR;U_OCR];
%     Ut_E=[Ut_E;U_E];
    
end

%disp(counter_rand);
if length_20
    disp('NOW STARTING WITH POLYMER >= THAN 20!!!')
    for i=t_found+1:t_found+t
        %disp(i);
        is_dsRNA=false;
          
            %% OBTAIN INITIAL IMPORTANT VALUES: Lsize, W, P1, P2, R1, Pn and Pi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Lsize=size(L,2); % Number of colums in L (i.e., the number of molecules that are in the pool at that time)
            W=zeros(1,Lsize); % Weigth of the molecules of the pool, that is proportional to their length.
    
            for ele=1:Lsize
                % At time step 1 they are all monomers, so they should all weight the same: 1 
                if any(L{1,ele}=='_') % If there is a dsRNA
                    subs=size(find(L{1,ele}=='-'),2);
                    W(ele)=size(L{1,ele},2)-1-subs; % We substract the '_' if it is a dsRNA and also the possible empty spaces ('-'). E.g. L{1,ele}='gcg_-g-', then W(ele) = 7-1-2 = 4 nts.
                else
                    W(ele)=size(L{1,ele},2);
                end
    
            end
    
            Wn=W./sum(W); % Normalization
        %     reset(RandStream.getGlobalStream,sum(100*clock)) 
            P1=sum(rand >= cumsum([0,Wn])); % Random selection of one molecule from the pool as a function of their weigths
            P2=randi(Lambda); % Randomly selected position where the nt or polymer is going to fall
    
    
            R1_L=L{1,P1}; % Extracting the Nt or polymer that is going to fall to the clay
    
            if any(R1_L=='_') % If we have a dsRNA, choose the strand that is going to fall to the clay => Always the one that was facing the clay (the molecule cannot turn around).
                sep=find(R1_L=='_');
                R1=R1_L(1:(sep-1)); % Left strand is always going to be the one that was in direct touch with clay
                R1_U=R1_L((sep+1):end); % Extract also the sequence of the complementary strand (right one).
    
                % Extract the OCR label for both strands
                OCR_both=L{2,P1};
                E_both=L{4,P1};
                OCR_left=OCR_both(1:(sep-1));
                E_left=E_both(1:(sep-1));            
                OCR_right=OCR_both((sep+1):end);
                E_right=E_both((sep+1):end);
                OCR=OCR_left;
                E=E_left;
    
    
                is_dsRNA=true;        
            else
                R1=R1_L;
                OCR=L{2,P1};
                E=L{4,P1};
            end
    
            % We choose randomly which of the nts of R1 molecule interacts with the clay:
            % Pn= position of R1 that interacts, Pi= position(s) of clay with which R1 interacts.
    
            if size(R1,2)==1 % The selected molecule is a single nt
            Pi=P2; % Place of interaction  
            Pn=1; 
            else % The selected molecule is an oligomer
            Pn=randi(size(R1,2)); % Pn is the position (nt) of the R1 chain that interacts with the position P2 
            Pi=(P2-Pn+1):(P2+size(R1,2)-Pn); % Complete place of interaction: it is a vector that indicates all the clay positions with which R1 interacts
            end
    
            neg=any(Pi<=0); % the position is out in the left edge
            % (P2+size(R1,2)-Pn) <= size(O_C,2) % <=== is the condition for the right edge
            
            
            %% EVALUATE STRUCTURE POSSIBLE FIXATION TO LEVEL 0 OR 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            if (P2+size(R1,2)-Pn) <= size(O_C,2) & neg==0 % Conditions in both edges: only if they fall inside the available space (Lambda) they can get adsorbed
    
                if sum(O_C(Pi))==0 % It falls in clay (level 0) with available place
    
                    Rup=rupture(alpha(i),size(R1,2)); % Evaluate desorption probability using the helper function rupture.m
                    counter_rand=counter_rand+1;
    
                    if Rup==0 % If links survive to one rupture round, actualize the variables
    
                        O_C(Pi)=ones(size(Pi));
                        C(Pi)=R1;
                        C_OCR(Pi)=OCR;
                        C_bl(Pi)=L{3,P1};
                        C_E(Pi)=E; 
    
                        if Pi(1)~=1 % Evaluate the polymerization of R1 with the nt to its left (if there is one)
                            if C_bl(Pi(1)-1)~='-'
                                if C_OCR(Pi(1))=='C' & C_OCR(Pi(1)-1)=='C'
                                    C_bl(Pi(1)-1)='1';
                                else
                                    C_bl(Pi(1)-1)='0';
                                end
                            end
                        end
    
                        if Pi(end)~=Lambda % Evaluate the polymerization of R1 with the nt to its right (if there is one)
                            if C_bl(Pi(end)+1)~='-'
                                if C_OCR(Pi(end))=='C' & C_OCR(Pi(end)+1)=='C'
                                    C_bl(Pi(end))='1';
                                else
                                    C_bl(Pi(end))='0';
                                end
                            end
                        end
    
                        if is_dsRNA % If it is a dsRNA we also have to update U, U_ORC and O_U in addition to O_C, C, C_OCR and C_bl
    
                            U(Pi)=R1_U;
                            Pi_U=Pi(R1_U~='-'); % Occupied place by R1_U, because it can be an oligo that is smaller than its opposite strand. E.g R1_L='gcg_-g-'. Here if Pi is for example [4 5 6], Pi_U would be just [5]
                            O_U(Pi_U)=ones(size(Pi_U));                
                            U_OCR(Pi)=OCR_right;
                            U_E(Pi)=E_right;
    
                        end            
    
                        % Delete the falling molecule from the pool   
                        L=reduction_E1(L,P1); 
                        total_fixed_0=total_fixed_0+1;
                        array_fixed_0(size(R1,2))=array_fixed_0(size(R1,2))+1;
    
                    end
    
                elseif sum(O_C(Pi))==size(R1,2) % If all positions at level 0 are occupied, but none at level 1 is occupied, R1 can fall in level 1 only if it is ssRNA.
    
                    if sum(O_U(Pi))==0 & is_dsRNA==false % Necessary conditions: available place at level 1 and only ssRNA can hybridize to level 1
    
                        % 1st) EVALUATE COMPLEMENTARITY
    
                        R2=C(Pi); % RNA molecule in clay against which we should evaluate the complementarity of R1
                        RC=complementary2(D,R1,R2); % Evaluate nucleobase complementarity according to rules set on D and using helper funtion complementary2.m
                        
                        if RC==1 % If there is full complementarity
                             
                             % 2nd) EVALUATE DENATURATION PROBABILITY
    
                             Rup = rupture(beta(i),size(R1,2));
                             counter_rand=counter_rand+1;
    
                             if Rup==0 % If links survive to one rupture round, actualize the variables
    
                                 U(Pi)=R1;
                                 O_U(Pi)=ones(size(Pi));
    
                                 % 3rd) CHANGE THE OCR AND E LABELS OF R1  ACCORDING TO RULES IN OCR_table_E1.m and E_table_E1_V2.m
    
                                 R1_OCR=L{2,P1}; % OCR labels of the molecule R1 that fell to template
                                 R1_E=L{4,P1};
                                 R2_OCR=C_OCR(Pi); % OCR labels of the template R2 above which R1 fell
                                 R2_E=C_E(Pi);
                                 new_R1_OCR=OCR_table_E1(R1_OCR,R2_OCR); 
                                 U_OCR(Pi)=new_R1_OCR; 
                                 new_R1_E=E_table_E1_V2(R1_E,R2_E);
                                 U_E(Pi)=new_R1_E;
    
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
    
    
            %% EVALUATE RUPTURE PROBABILITIES OF ALL THE MOLECULES IN LEVEL 0, 1 OR 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % -- The order is not important, so we are going to start testing level 1 and 2 (pool) linkages (hydrogen bonds)
            
            % LEVEL 2 (POOL) => Evaluate the denaturation probability of all dsRNA molecules present in the pool in the current time step
            for pos=1:size(L,2)
                
                if any(L{1,pos}=='_') % If we have dsRNA, evaluate denaturation probability
            
                    ds_oligo=L{1,pos};
                    ds_OCR=L{2,pos};
                    ds_E=L{4,pos};
                    bl_l=L{3,pos}; % Cause we only store the bl label of level 0
                    separ=find(ds_oligo=='_');
                    oligo_r=ds_oligo((separ+1):end); % level 1 oligomer (to the right of '_' separator)
                    OCR_r=ds_OCR((separ+1):end);  
                    E_r=ds_E((separ+1):end); 
                    oligo_l=ds_oligo(1:(separ-1)); % level 0 oligomer (to the left of '_' separator)
                    OCR_l=ds_OCR(1:(separ-1));
                    E_l=ds_E(1:(separ-1));
    
                    empty_pos=find(oligo_r=='-');
                    empty_pos=[0,empty_pos,length(oligo_r)+1];
                    diff_empty_pos=diff(empty_pos)-1;
    %                 diff_empty_pos=diff(empty_pos);
                    reduc=false; 
    
                    if size(find(diff_empty_pos>0),2)>1 % If we have for example: oligo_r = '---auuagg----cgaguag-'
                                            
                        % --- 1st) Add new entry for the right ssRNA products that are the products (in plural) of dsRNA rupture
                        
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
                                            nE_r=E_r(1:k(1)-1);
                                            if any(nOCR_r=='R')
                                                R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r};                        
                                            end
    
                                            if any(OCR_r=='C')
                                                C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r};
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
                                            if k(1)-1==1
                                                nE_r='!';
                                            end
                                            L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                            oligo_r(1:k(1)-1)='-';
                                            OCR_r(1:k(1)-1)='-';
                                            E_r(1:k(1)-1)='-';
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
                                            nE_r=E_r((k(z)+1):(k(z+1)-1));
                                            
                                            if any(nOCR_r=='R')
                                                R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r};                        
                                            end
    
                                            if any(OCR_r=='C')
                                                C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r};
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
                                            if (k(z+1)-k(z)-1)==1
                                                nE_r='!';
                                            end
                                            L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                            oligo_r((k(z)+1):(k(z+1)-1))='-';
                                            OCR_r((k(z)+1):(k(z+1)-1))='-';
                                            E_r((k(z)+1):(k(z+1)-1))='-';
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
                                            nE_r=E_r(k(temps)+1:size(oligo_l,2));
    
                                            if any(nOCR_r=='R')
                                                R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r};                        
                                            end
    
                                            if any(OCR_r=='C')
                                                C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r};
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
                                            if (size(oligo_l,2)-k(temps))==1
                                                nE_r='!';
                                            end
                                            L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                            oligo_r(k(temps)+1:size(oligo_l,2))='-';
                                            OCR_r(k(temps)+1:size(oligo_l,2))='-';
                                            E_r(k(temps)+1:size(oligo_l,2))='-';
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
                                            nE_r=E_r((k(z)+1):(k(z+1)-1));
    
                                            if any(nOCR_r=='R')
                                                R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r};                        
                                            end
    
                                            if any(OCR_r=='C')
                                                C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r};
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
                                            if (k(z+1)-k(z)-1)==1
                                                nE_r='!';
                                            end
                                            L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                            oligo_r((k(z)+1):(k(z+1)-1))='-';
                                            OCR_r((k(z)+1):(k(z+1)-1))='-';
                                            E_r((k(z)+1):(k(z+1)-1))='-';
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
                                    nE_r=E_r(1:k(1)-1);
    
                                    if any(nOCR_r=='R')
                                        R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r};                        
                                    end
    
                                    if any(OCR_r=='C')
                                        C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r};
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
                                    if (k(1)-1)==1
                                        nE_r='!';
                                    end
                                    L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                    oligo_r(1:k(1)-1)='-';
                                    OCR_r(1:k(1)-1)='-';
                                    E_r(1:k(1)-1)='-';
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
                                    nE_r=E_r(k(temps)+1:size(oligo_l,2));
    
                                    if any(nOCR_r=='R')
                                        R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r};                        
                                    end
    
                                    if any(OCR_r=='C')
                                        C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r};
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
                                    if (size(oligo_l,2)-k(temps))==1
                                        nE_r='!';
                                    end
                                    L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                    oligo_r(k(temps)+1:size(oligo_l,2))='-';
                                    OCR_r(k(temps)+1:size(oligo_l,2))='-';
                                    E_r(k(temps)+1:size(oligo_l,2))='-';
                                    total_stripped_2=total_stripped_2+1;
                                    array_stripped_2(size(oligo_l,2)-k(temps))=array_stripped_2(size(oligo_l,2)-k(temps))+1;
                                    reduc=true;
                                end
                            end
    
    
                        else
                             %(if temps==0), it is not necessary to evaluate the entire chain rupture here, because we already evaluated it
                        end
    
                        if reduc
                            
                            % --- 2nd) Add new entry for the RNA that is left, that can be just the left ssRNA or another dsRNA
                            if all(oligo_r=='-')
                                if size(E_l,2)==1
                                    E_l='!';
                                end
                                L=ampliation_E1_MODIF(L,oligo_l,OCR_l,bl_l,E_l);
                            else
                                L=ampliation_E1_MODIF(L,strcat(oligo_l,'_',oligo_r),strcat(OCR_l,'_',OCR_r),bl_l,strcat(E_l,'_',E_r));
                            end
    
                            % --- 3rd) Delete the entry in L cell of the original dsRNA
                            isdsOligo = cellfun(@(x)isequal(x,ds_oligo),L(1,:));
                            [row,col] = find(isdsOligo); % Find the original dsRNA position, because we have added new entries to L, and positions changed
                            L=reduction_E1(L,col); 
    
                        end
    
                    else % If oligo_r is simpler
    
                        size_rup=size(find(oligo_r~='-'),2); 
                        
                        if size_rup~=0
                            Rup=rupture(beta(i),size_rup);
                            counter_rand=counter_rand+1;
                            
                            if Rup==1 
                                noligo_r=oligo_r(oligo_r~='-'); % Delete possible empty ('-') positions of right oligo
                                nOCR_r=OCR_r(OCR_r~='-');
                                nE_r=E_r(E_r~='-');
                                if size(nE_r,2)==1
                                    nE_r='!';
                                end
    
                                if any(OCR_r=='R')
                                    R_CELL(end+1,:)={OCR_r, oligo_r, bl_l, E_r};                        
                                end
                                
                                if any(OCR_r=='C')
                                    C_CELL(end+1,:)={OCR_r, oligo_r, bl_l, E_r};
                                end
                                
                                % -- Actualize L (the pool)
                                % --- 1st) Delete the entry in L cell of the dsRNA
                                L=reduction_E1(L,pos); 
                                % --- 2nd) Add two new entries for the two ssRNA that are the products dsRNA rupture
                                if size(E_l,2)==1
                                    E_l='!';
                                end
                                L=ampliation_E1_MODIF(L,oligo_l,OCR_l,bl_l,E_l);
                                bl_r='';
                                if size(noligo_r,2)==1
                                    bl_r='X';
                                else
                                    bl_r(1:size(noligo_r,2))='0';
                                    bl_r(end)='X';
                                end
                                new_OCR_r=nOCR_r;
                                new_OCR_r(nOCR_r=='R')='O';
                                L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r);
                                
                            end
                        end
                    end
                    
                end
                
            end
            
            % LEVEL 1 => Evaluate the denaturation probability of the molecules hybridized to level 1 in the current time step
            
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
                                E_chain=U_E(1:k(1)-1);
    
                                if any(OCR_chain=='R')
                                    R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};                        
                                end
    
                                if any(OCR_chain=='C')
                                    C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};
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
                                if (k(1)-1)==1
                                    E_chain='!';
                                end
                                L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                                U(1:k(1)-1)='-';
                                O_U(1:k(1)-1)=0;
                                U_OCR(1:k(1)-1)='-';
                                U_E(1:k(1)-1)='-';
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
                                E_chain=U_E((k(z)+1):(k(z+1)-1));
    
                                if any(OCR_chain=='R')
                                    R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten,E_chain};                        
                                end
    
                                if any(OCR_chain=='C')
                                    C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten,E_chain};
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
    
                                if (k(z+1)-k(z)-1)==1
                                    E_chain='!';
                                end
                                L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                                U((k(z)+1):(k(z+1)-1))='-';
                                O_U((k(z)+1):(k(z+1)-1))=0;
                                U_OCR((k(z)+1):(k(z+1)-1))='-';
                                U_E((k(z)+1):(k(z+1)-1))='-';
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
                                E_chain=U_E(k(temps)+1:Lambda);
    
                                if any(OCR_chain=='R')
                                    R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};                        
                                end
    
                                if any(OCR_chain=='C')
                                    C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};
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
    
                                if (Lambda-k(temps))==1
                                    E_chain='!';
                                end
                                L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                                U(k(temps)+1:Lambda)='-';
                                O_U(k(temps)+1:Lambda)=0;
                                U_OCR(k(temps)+1:Lambda)='-';
                                U_E(k(temps)+1:Lambda)='-';
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
                                E_chain=U_E((k(z)+1):(k(z+1)-1));                           
    
                                if any(OCR_chain=='R')
                                    R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};                        
                                end
    
                                if any(OCR_chain=='C')
                                    C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};
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
    
                                if (k(z+1)-k(z)-1)==1
                                    E_chain='!';
                                end
                                L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                                U((k(z)+1):(k(z+1)-1))='-';
                                O_U((k(z)+1):(k(z+1)-1))=0;
                                U_OCR((k(z)+1):(k(z+1)-1))='-';
                                U_E((k(z)+1):(k(z+1)-1))='-';
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
                        E_chain=U_E(1:k(1)-1);
    
                        if any(OCR_chain=='R')
                            R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};                        
                        end
    
                        if any(OCR_chain=='C')
                            C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain};
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
                        if (k(1)-1)==1
                            E_chain='!';
                        end
                        L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                        U(1:k(1)-1)='-';
                        O_U(1:k(1)-1)=0;
                        U_OCR(1:k(1)-1)='-';
                        U_E(1:k(1)-1)='-';
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
                        E_chain=U_E(k(temps)+1:Lambda);
                        if any(OCR_chain=='R')
                            R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};                        
                        end
    
                        if any(OCR_chain=='C')
                            C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};
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
    
                        if (Lambda-k(temps))==1
                            E_chain='!';
                        end
                        L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                        U(k(temps)+1:Lambda)='-';
                        O_U(k(temps)+1:Lambda)=0;
                        U_OCR(k(temps)+1:Lambda)='-';
                        U_E(k(temps)+1:Lambda)='-';
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
                    E_chain=U_E;
                    if any(OCR_chain=='R')
                        R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};                        
                    end
    
                    if any(OCR_chain=='C')
                        C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain};
                    end
    
                    bl_new='';
                    bl_new(1:size(nt_chain,2))='0'; % All links are "good", because they were generated by template depending polymerization
                    bl_new(end)='X';
                    new_OCR_chain=OCR_chain;
                    new_OCR_chain(OCR_chain=='R')='O';   
    
                    L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                    U(1:Lambda)='-';
                    O_U(1:Lambda)=0;
                    U_OCR(1:Lambda)='-';
                    U_E(1:Lambda)='-';
                    total_stripped_1=total_stripped_1+1;
                    array_stripped_1(size(O_U,2))=array_stripped_1(size(O_U,2))+1;
                end
                
            end
    
            % LEVEL 0 => Evaluate the desorption probability of all the molecules adsorbed to level 0
    
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
                                    E_chain=strcat(C_E(1:k(1)-1),'_',U_E(1:k(1)-1));
                                    U(1:k(1)-1)='-';
                                    O_U(1:k(1)-1)=0;
                                    U_OCR(1:k(1)-1)='-';
                                    U_E(1:k(1)-1)='-';
                                else
                                    OCR_chain=C_OCR(1:k(1)-1);
                                    nt_chain=C(1:k(1)-1);
                                    bl_chain=C_bl(1:k(1)-1);  
                                    
                                    if (k(1)-1)==1
                                        E_chain='!';
                                    else
                                        E_chain=C_E(1:k(1)-1);
                                    end
                                end
                         
                                
                                %-- Actualize the variables
                                L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                                C(1:k(1)-1)='-';
                                O_C(1:k(1)-1)=0;
                                C_OCR(1:k(1)-1)='-';
                                C_bl(1:k(1)-1)='-';
                                C_E(1:k(1)-1)='-';
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
                                    E_chain=strcat(C_E((k(z)+1):(k(z+1)-1)),'_',U_E((k(z)+1):(k(z+1)-1)));
                                    U((k(z)+1):(k(z+1)-1))='-';
                                    O_U((k(z)+1):(k(z+1)-1))=0;
                                    U_OCR((k(z)+1):(k(z+1)-1))='-';
                                    U_E((k(z)+1):(k(z+1)-1))='-';
                                else
                                    OCR_chain=C_OCR((k(z)+1):(k(z+1)-1));
                                    nt_chain=C((k(z)+1):(k(z+1)-1));
                                    bl_chain=C_bl((k(z)+1):(k(z+1)-1)); 
                                     
                                    if (k(z+1)-k(z)-1)==1
                                        E_chain='!';
                                    else
                                        E_chain=C_E((k(z)+1):(k(z+1)-1));
                                    end
                                end
                         
                                
                                %-- Actualize the variables
                                L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                                C((k(z)+1):(k(z+1)-1))='-';
                                O_C((k(z)+1):(k(z+1)-1))=0;
                                C_OCR((k(z)+1):(k(z+1)-1))='-';
                                C_bl((k(z)+1):(k(z+1)-1))='-';
                                C_E((k(z)+1):(k(z+1)-1))='-';
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
                                    E_chain=strcat(C_E(k(temps)+1:Lambda),'_',U_E(k(temps)+1:Lambda));
                                    U(k(temps)+1:Lambda)='-';
                                    O_U(k(temps)+1:Lambda)=0;
                                    U_OCR(k(temps)+1:Lambda)='-';
                                    U_E(k(temps)+1:Lambda)='-';
                                else
                                    OCR_chain=C_OCR(k(temps)+1:Lambda);
                                    nt_chain=C(k(temps)+1:Lambda);
                                    bl_chain=C_bl(k(temps)+1:Lambda);  
                                    if (Lambda-k(temps))==1
                                        E_chain='!';
                                    else
                                        E_chain=C_E(k(temps)+1:Lambda);  
                                    end
                                end
                         
                                
                                %-- Actualize the variables
                                L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                                C(k(temps)+1:Lambda)='-';
                                O_C(k(temps)+1:Lambda)=0;
                                C_OCR(k(temps)+1:Lambda)='-';
                                C_bl(k(temps)+1:Lambda)='-';
                                U_E(k(temps)+1:Lambda)='-';
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
                                    E_chain=strcat(C_E((k(z)+1):(k(z+1)-1)),'_',U_E((k(z)+1):(k(z+1)-1)));
                                    U((k(z)+1):(k(z+1)-1))='-';
                                    O_U((k(z)+1):(k(z+1)-1))=0;
                                    U_OCR((k(z)+1):(k(z+1)-1))='-';
                                    U_E((k(z)+1):(k(z+1)-1))='-';
                                else
                                    OCR_chain=C_OCR((k(z)+1):(k(z+1)-1));
                                    nt_chain=C((k(z)+1):(k(z+1)-1));
                                    bl_chain=C_bl((k(z)+1):(k(z+1)-1)); 
                                    if (k(z+1)-k(z)-1)==1
                                        E_chain='!';
                                    else
                                        E_chain=C_E((k(z)+1):(k(z+1)-1));
                                    end
                                end
                         
                                
                                %-- Actualize the variables
                                L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                                C((k(z)+1):(k(z+1)-1))='-';
                                O_C((k(z)+1):(k(z+1)-1))=0;
                                C_OCR((k(z)+1):(k(z+1)-1))='-';
                                C_bl((k(z)+1):(k(z+1)-1))='-';
                                C_E((k(z)+1):(k(z+1)-1))='-';
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
                            E_chain=strcat(C_E(1:k(1)-1),'_',U_E(1:k(1)-1));
                            U(1:k(1)-1)='-';
                            O_U(1:k(1)-1)=0;
                            U_OCR(1:k(1)-1)='-';
                            U_E(1:k(1)-1)='-';
                        else
                            OCR_chain=C_OCR(1:k(1)-1);
                            nt_chain=C(1:k(1)-1);
                            bl_chain=C_bl(1:k(1)-1); 
                            if (k(1)-1)==1
                                E_chain='!';
                            else
                                E_chain=C_E(1:k(1)-1);
                            end
                        end
                 
                        
                        %-- Actualize the variables
                        L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                        C(1:k(1)-1)='-';
                        O_C(1:k(1)-1)=0;
                        C_OCR(1:k(1)-1)='-';
                        C_bl(1:k(1)-1)='-';
                        C_E(1:k(1)-1)='-';
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
                            E_chain=strcat(C_E(k(temps)+1:Lambda),'_',U_E(k(temps)+1:Lambda));
                            U(k(temps)+1:Lambda)='-';
                            O_U(k(temps)+1:Lambda)=0;
                            U_OCR(k(temps)+1:Lambda)='-';
                            U_E(k(temps)+1:Lambda)='-';
                        else
                            OCR_chain=C_OCR(k(temps)+1:Lambda);
                            nt_chain=C(k(temps)+1:Lambda);
                            bl_chain=C_bl(k(temps)+1:Lambda);  
                            if (Lambda-k(temps))==1
                                E_chain='!';
                            else
                                E_chain=C_E(k(temps)+1:Lambda);
                            end
                        end
                 
                        
                        %-- Actualize the variables
                        L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                        C(k(temps)+1:Lambda)='-';
                        O_C(k(temps)+1:Lambda)=0;
                        C_OCR(k(temps)+1:Lambda)='-';
                        C_bl(k(temps)+1:Lambda)='-';
                        C_E(k(temps)+1:Lambda)='-';
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
                        E_chain=strcat(C_E,'_',U_E);
                        U(1:Lambda)='-';
                        O_U(1:Lambda)=0;
                        U_OCR(1:Lambda)='-';
                        U_E(1:Lambda)='-';
                    else
                        OCR_chain=C_OCR;
                        nt_chain=C;
                        bl_chain=C_bl; 
                        E_chain=C_E;
                    end
             
                    
                    %-- Actualize the variables
                    L=ampliation_E1_MODIF(L,nt_chain,OCR_chain,bl_chain,E_chain);
                    C(1:Lambda)='-';
                    O_C(1:Lambda)=0;
                    C_OCR(1:Lambda)='-';
                    C_bl(1:Lambda)='-';
                    C_E(1:Lambda)='-';
                    total_stripped_0=total_stripped_0+1;
                    array_stripped_0(size(O_C,2))=array_stripped_0(size(O_C,2))+1;
                end
            end
    
            if i==t
                snapshot_tmax={};
                % Calculate activity cell to store it:
                        Activity_tmax=cell(5,3);
                        Activity_tmax{1,1}='fixed level 0';
                        Activity_tmax{2,1}='stripped level 0';
                        Activity_tmax{1,2}=total_fixed_0;
                        Activity_tmax{2,2}=total_stripped_0;
                        Activity_tmax{1,3}=array_fixed_0;
                        Activity_tmax{2,3}=array_stripped_0;
                        Activity_tmax{3,1}='fixed level 1';
                        Activity_tmax{4,1}='stripped level 1';
                        Activity_tmax{3,2}=total_fixed_1;
                        Activity_tmax{4,2}=total_stripped_1;
                        Activity_tmax{3,3}=array_fixed_1;
                        Activity_tmax{4,3}=array_stripped_1;
                        Activity_tmax{5,1}='stripped level 2';
                        Activity_tmax{5,2}=total_stripped_2;
                        Activity_tmax{5,3}=array_stripped_2;
               snapshot_tmax{1}=Activity_tmax;
               snapshot_tmax{2}=C_CELL;
               snapshot_tmax{3}=R_CELL;
               snapshot_tmax{4}=C;
               snapshot_tmax{5}=C_OCR;
               snapshot_tmax{6}=C_bl;
               snapshot_tmax{7}=C_E;
            end
    
    % ---- Update variables. NOTICE: Uncomment following lines to obtain all possible data from simulations in addition to the Activity, R_CELL and C_CELL variables (althoug notice that this will require additional time per simulation)----
    
    %         % Update the pool and the occupancy arrays
    %         Ltemp='|';
    %         LOCRtemp='|';
    %         Lbltemp='|';
    %         LEtemp='|';
    %         for n=1:size(L,2)
    %             tempL=[L{1,n},'|'];
    %             tempLOCR=[L{2,n},'|'];
    %             tempLbl=[L{3,n},'|'];
    %             tempLE=[L{4,n},'|'];
    %             tempL2=[Ltemp,tempL];
    %             tempLOCR2=[LOCRtemp,tempLOCR];
    %             tempLbl2=[Lbltemp,tempLbl];
    %             tempLE2=[LEtemp,tempLE];
    %             Ltemp=tempL2; %Ltemp has a structure like this: '|A|A|T|AG|A|AAC|...', where | separates the different oligos/nts.
    %             LOCRtemp=tempLOCR2; %LOCRtemp has a structure like this: '|O|O|O|CC|O|OCC|...', where | separates the different oligos/nts.
    %             Lbltemp=tempLbl2; %Lbltemp has a structure like this: '|0|0|0|0X|0|00X|...', where | separates the different oligos/nts.
    %             LEtemp=tempLE2; %Lbltemp has a structure like this: '|!|!|!|AB|!|!JK|...', where | separates the different oligos/nts.
    %         end
    % 
             tf=i; % Utilitarian variable
    % 
    %     if i~=1
    % 
    %         %% -- Lt
    % 
    %         if size(Ltemp,2)==size(Lt(end,:),2)
    %             Lt=[Lt;Ltemp];
    %         else
    %             prev_size=size(Lt(end,:),2);
    %             new_size=size(Ltemp,2);
    % 
    %             if prev_size>new_size
    %                 Ltemp(end+1:size(Lt(end,:),2))='|';
    %             elseif new_size>prev_size
    %                 Lt(1:end,end+1:size(Ltemp(end,:),2))='|';
    %             end
    %             Lt=[Lt;Ltemp];
    %         end
    % 
    %          %% -- LtOCR
    %         
    %         if size(LOCRtemp,2)==size(Lt_OCR(end,:),2)
    %             Lt_OCR=[Lt_OCR;LOCRtemp];
    %         else
    %             prev_size=size(Lt_OCR(end,:),2);
    %             new_size=size(LOCRtemp,2);
    % 
    %             if prev_size>new_size
    %                 LOCRtemp(end+1:size(Lt_OCR(end,:),2))='|';
    %             elseif new_size>prev_size
    %                 Lt_OCR(1:end,end+1:size(LOCRtemp(end,:),2))='|';
    %             end
    %             Lt_OCR=[Lt_OCR;LOCRtemp];
    %         end
    % 
    %         %% -- Ltbl
    % 
    %         if size(Lbltemp,2)==size(Lt_bl(end,:),2)
    %             Lt_bl=[Lt_bl;Lbltemp];
    %         else
    %             prev_size=size(Lt_bl(end,:),2);
    %             new_size=size(Lbltemp,2);
    % 
    %             if prev_size>new_size
    %                 Lbltemp(end+1:size(Lt_bl(end,:),2))='|';
    %             elseif new_size>prev_size
    %                 Lt_bl(1:end,end+1:size(Lbltemp(end,:),2))='|';
    %             end
    %             Lt_bl=[Lt_bl;Lbltemp];
    %         end
    % 
    %         %% -- LtE
    % 
    %         if size(LEtemp,2)==size(Lt_E(end,:),2)
    %             Lt_E=[Lt_E;LEtemp];
    %         else
    %             prev_size=size(Lt_E(end,:),2);
    %             new_size=size(LEtemp,2);
    % 
    %             if prev_size>new_size
    %                 LEtemp(end+1:size(Lt_E(end,:),2))='|';
    %             elseif new_size>prev_size
    %                 Lt_E(1:end,end+1:size(LEtemp(end,:),2))='|';
    %             end
    %             Lt_E=[Lt_E;LEtemp];
    %         end
    % 
    % 
    %     else
    %         Lt=[Lt;Ltemp];
    %         Lt_OCR=[Lt_OCR;LOCRtemp];
    %         Lt_bl=[Lt_bl;Lbltemp];
    %         Lt_E=[Lt_E;LEtemp];
    % 
    %     end
    
%     Ot_C=[Ot_C;O_C];
%     Ot_U=[Ot_U;O_U];
%     Ct=[Ct;C];
%     Ct_OCR=[Ct_OCR;C_OCR];
    % Ct_bl=[Ct_bl;C_bl];
%     Ct_E=[Ct_E;C_E];
%     Ut=[Ut;U];
%     Ut_OCR=[Ut_OCR;U_OCR];
%     Ut_E=[Ut_E;U_E];
    
    
    end
    
end

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