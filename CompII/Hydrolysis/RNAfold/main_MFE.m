function [O_seq,Lt,Lt_OCR,Lt_bl,Lt_E,tf,Ut,Ut_OCR,Ut_E,Ct,Ct_OCR,Ct_bl,Ct_E,Ot_C,Ot_U,Activity,R_CELL,C_CELL, HYDROLYZED, N_pbonds, H_pbonds, folded_how_many,mfe_how_many, total_pols_pool, mean_energy, total_denat_fold, total_pols_C,total_pols_U]=main_MFE(A,D,N,positions,t,beta,alpha,polymer_size,k_hyd, seq) 
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

Lambda=positions; % Number of simulated clay positions (usually Lambda=200)

L=cell(5,sum(N)); % Create a cell to store the composition of the pool, that is a matrix like structure than can harbor any data type. Here they are all strings. 
                  % sum(N) is the total amount of nts in the pool at the beginning of the simulation. 
                  % Each column contains one molecule and its related labels.
                  % Row 1 => The nt sequence of the ssRNAs or dsRNAs (dsRNA chains are separated by '_'). E.g {a} {u} {auu_ua} {cgcg} <= Example of 4 molecules present in the pool, two are
                  % single nucleotides, the third one is a dsRNA with 3 nt in one strand and 2 nt in the complementary strand, the fourth is a 4-nt long ssRNA.
                  % Row 2 => The OCR label of each nt of the molecule. E.g {O} {O} {OOO_CC} {CCCC} % Notice that O here is a letter, not the number
                  % Row 3 => The bl label of the molecule. E.g {X} {X} {00X_0X} {010X} % The X means that there isn't a link (because the number of links in a molecule is their length - 1),
                  % also, the sequence starts at the first link, but it has the same length as the sequence of nts, therefore, last element will always be an X.
                  % Row 4 => The E label of each nt of the molecule. E.g {!} {!} {AB!_AB} {!!JK}
                  % Row 5 => MFE, folding energy and folding state. It is a 1x3 cell that stores the folding in bracket notation, the energy and the folding state in that time step. E.g.{{'..(((..)))..',[-0.7],[1]}} {{'.....',[0],[1]}}. If
                  % the folding state is [1] it means that it is folded in that time step and if it is [0] it means that it has denaturalized in the previous time step when evaluating Prup. 
                  % It should be refolded again each time step and reevaluated its Prup. The RNAs whose MFE is already the unfolded strand, always will have parameter [0] here for convenience.


%Lt=''; % Created to store the info of row 1 of L variable (pool) for each time step.
% Lt_OCR=''; % Created to store the info of row 2 of L variable for each time step.
% Lt_bl=''; % Created to store the info of row 3 of L variable for each time step.
% Lt_E=''; % Created to store the info of row 4 of L variable for each time step.
Lt=cell(t,1); 
Lt_OCR=cell(t,1);
Lt_bl=cell(t,1);
Lt_E=cell(t,1);
% Lt_fold=cell(t,1);

%folded_how_many=NaN(t,1);
%total_pols_pool=Nan(t,1);
%mfe_how_many=Nan(t,1);
folded_how_many=cell(t,2); % {total folded, [their lengths]}
mfe_how_many=cell(t,2); % {total with mfe, [their lengths]}.
total_pols_pool=cell(t,2); % {total pols, [their lengths]}.
mean_energy=NaN(t,1);
total_denat_fold=0;
total_pols_C=cell(t,2);
total_pols_U=cell(t,2); 


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
for i=1:size(N,2) % For each type of nucleotide (4 normally)
    for n=1:N(i) % For each of the nucleotides of the same kind (2*template lenght normally)
        L{1,l}=A(i);
        L{2,l}='O';
        L{3,l}='X';
        L{4,l}='!'; 
        L{5,l}={'.',[],[0]};
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

% -- Place original polymer (sequence) in the pool

ascii='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';

if isempty(seq) % -- If not predetermined sequence is provided, the original sequence is random
    nums = randi(size(A,2),[1 polymer_size]);
    seq = A(nums);
end

O_seq=seq;

L{1,end+1}=seq;
L{2,end}=repmat('O',[1,polymer_size]);
L{3,end}=[repmat('0',[1,polymer_size-1]),'X'];
L{4,end}=ascii(1:polymer_size);
cmd = ['echo ', seq, ' | RNAfold'];
[status, result] = system(cmd);
clean_result = regexprep(result, '\x1B\[[0-9;]*[A-Za-z]', '');
lines = strsplit(strtrim(clean_result), '\n');
struct_energy=split(lines{2});
struct=struct_energy{1};
energy=struct_energy{3};
energy=str2num(energy(1:end-1));
if all(struct=='.')
    L{5,end}={struct,[energy],[0]};
else
    L{5,end}={struct,[energy],[1]};
end

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
R_CELL={}; % Every row stores a molecule with at least one nt being of type 'R' and all of its associated labels (OCR, bl and E), and its folding
C_CELL={}; % Every row stores a molecule with at least one nt being of type 'C' and all of its associated labels (OCR, bl and E), and its folding

% -- Variable to store the molecules that hydrolyze, the new molecules formed by their hydrolysis and the time step in which this occurs
HYDROLYZED={}; % Every row stores a molecule, and the columns mean: 
               % "broken molecule sequence" "OCR_label" "bl_label" "E_label" {cell array with molecules formed after hydrolysis} "time step"
% -- Variable to store the total number of susceptible p-bonds at each time step
N_pbonds=[];
% -- Variable to store the total number of hydrolyzed p-bonds at each time step
H_pbonds=[];

counter_rand=0; % Utilitarian variable

%% MAIN CODE: Run a simulation over time

for i=1:t
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
            folded=0;
        else
            R1=R1_L;
            OCR=L{2,P1};
            E=L{4,P1};
            folding=L{5,P1};
            folded=folding{3};
        end

        if ~folded % Only the RNAs that are not folded can fall to clay.

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
        end

        if k_hyd(i)~=0
            %% EVALUATE HYDROLYSIS OF ssRNA IN THE POOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            p_bonds=cell(2,size(L,2)); % Cell to store the information of the covalent bonds (p-bonds) present in each molecule that are susceptible to hydrolysis.
                                       % Every new p-bond is going to be identified by a different integer number (first row) and we are going also to store the position of the p-bond in the molecule (second row).
                                       % Therefore if we have this three molecules in the pool: {a} {auu} {agg_uc-}, the p_bonds variable should be: {[0]} {[1 2]} {[3]} * integer id
                                       %                                                                                                             {[0]} {[1 2]} {[2]} * position of susceptible p-bond in molecule 
                                       %                                                                                                                                   (left strand for dsRNA molecules).
            
            pb_counter=1; % Counter of the total number of susceptible p-bonds in the pool at each time step.
    
            % 1st. CALCULATE THE P-BONDS THAT ARE SUSCEPTIBLE TO HYDROLYSIS IN EVERY MOLECULE AND ENUMERATE THEM
    
            for m=1:size(L,2) % Iterate through every molecule in the pool
    
                if any(L{1,m}=='_') % If there is a dsRNA
                    sep2=find(L{1,m}=='_');
                    left=L{1,m}(1:(sep2-1)); % Left strand is always going to be the one that was in direct touch with clay.
                    right=L{1,m}((sep2+1):end); % Extract also the sequence of the complementary strand (right one).
                    pb_array=[]; % Array to store the susceptible p-bonds ids of the molecule.
                    pb_m_array=[]; % Array to store the susceptible p-bonds positions in the molecule.
                    
                    % Notice that for the molecules in the pool that contain parts of dsRNA only the p-bonds that
                    % are not protected by two nucleotides forming Hydrogen bonds on top will be susceptible to hydrolysis:
                    
                    for nt=1:(size(left,2)-1) % Iterate trough each left molecule position to see the susceptible p-bonds based on the protection of complementary nucleotides in the right molecule
                        if right(nt)~='-' & right(nt+1)~='-' 
                            pb_array=[pb_array, 0]; % 0 means there is not a susceptible p-bond
                            pb_m_array=[pb_m_array, 0]; % 0 indicates here the position of a non-susceptible p-bond. Example of pb_m_array like this could be: [1 0 3], for molecule aaug_-ua-
                        else
                            pb_array=[pb_array, pb_counter]; % Include the integer of the new susceptible p-bond
                            pb_counter=pb_counter+1;
                            pb_m_array=[pb_m_array, nt]; % Position in molecule of the susceptible p-bond
                        end
                    end
                    p_bonds{1,m}=pb_array; % Notice pb_array and pb_m_array should have same length.
                    p_bonds{2,m}=pb_m_array;
    
                else % If the molecule is entirely a ssRNA
                    
                    m_folding=L{5,m};
                    m_folded=m_folding{3};
                    if m_folded % If molecule is folded, it cannot hydrolyze.
                        m_pb=0;    
                    else
                        m_size=size(L{1,m},2);
                        m_pb=m_size-1; % Number of susceptible p-bonds

                        % Return the molecule to the folded state if it has it:
                        m_struct=m_folding{1};
                        if any(m_struct~='.')
                            new_m_folding={m_folding{1},m_folding{2},[1]};
                            L{5,m}=new_m_folding;
                        end
                    end

                    if m_pb==0
                        p_bonds{1,m}=[0];
                        p_bonds{2,m}=[0];
                    else
                        new_pb_counter=pb_counter+m_pb;
                        p_bonds{1,m}=pb_counter:(new_pb_counter-1); % Generate new integer ids for all the p-bonds.
                        pb_counter=new_pb_counter;
                        p_bonds{2,m}=1:m_pb; % All p-bonds are susceptible in ssRNAs.
                    end
                end
    
            end
    
            % 2nd. CALCULATE THE NUMBER OF P-BONDS THAT ARE GOING TO BREAK IN THIS TIME STEP USING A BINOMIAL OF PARAMETERS N = total number of susceptible p-bonds, p = hydrolysis probability.
    
            num_breaks=binornd(pb_counter-1,k_hyd(i)); % Number of p-bonds that break
            
            N_pbonds=[N_pbonds, pb_counter-1];
            H_pbonds=[H_pbonds, num_breaks];
    
            if num_breaks ~=0
    
                % 3rd. RANDOMLY SELECT THE P-BONDS THAT ARE GOING TO BREAK AND ACTUALIZE THE L (POOL) VARIABLE

                %disp(k_hyd);
                %disp(strcat('pb_counter: ', num2str(pb_counter)));
                %disp(strcat('num_breaks: ', num2str(num_breaks)));
    
                pb_breaks=transpose(randsample(pb_counter-1,num_breaks));   % Specific p_bonds that break (integer identifiers)
                                                                            % randsample requires the Statistics toolbox, another function to be used istead could be randperm.
                                                                            % Notice we cannot use randi because we don't want to have repeated numbers.
    
                pb_breaks=sort(pb_breaks);
    
                % Find the specific p-bonds that are going to break and construct the new molecules after hydrolysis: 

                M_to_be_added={}; % Create cell to store the new molecules formed after hydrolysis, which are going to be added to the pool (L variable) at the end.        
    
                c={}; % Create empty cell in which we are going to link the p-bonds that will break to the molecules they belong to.
                      % c is a Nx2 cell in which each row stands for a molecule, in column 1 the molecule number as its identifier (idx, as it is sorted in L variable), 
                      % and in column 2 the position of the p-bonds of the molecule that will break (because we can have more than one p-bond hydrolysis in one molecule).
    
                for pb=pb_breaks 
    
                    idx = find(cellfun(@(x) any(x == pb), p_bonds(1,:))); % Index that corresponds to specific p-bond in p_bonds variable and to corresponding molecule in L variable
                    pb_ids = p_bonds{1,idx}; % Extract the entire array of p-bond integer identifiers for molecule idx
                    pb_pos = p_bonds{2,idx}; % Extract the entire array of p-bond positions in molecule idx
                    br = pb_pos(pb_ids==pb); % Extract the exact position of the p-bond that is going to break
    
                    if isempty(c)
                        c={idx,[br]};
                    else
                        c_key_idx = find(cellfun(@(x) any(x == idx), c(:,1))); % Index of row in c that corresponds to molecule number idx
    
                        if isempty(c_key_idx)
                            c(end+1,:)={idx, [br]};
                        else
                            c(c_key_idx,2)={[[c{c_key_idx,2}], br]}; 
                        end
                    end
    
                end
    
                M_to_be_removed=transpose(cell2mat(c(:,1))); % Array with the indexes of the molecules in the L variable that break, so we can remove them at the end all at once.
                for mol=1:size(c,1) % For every molecule that breaks, calculate the new molecules that are formed
                    molecule_idx=c{mol,1};
                    m_break = L{1, molecule_idx}; % Molecule to be broken
                    m_break_OCR = L{2,molecule_idx};
                    m_break_bl = L{3,molecule_idx};
                    m_break_E = L{4,molecule_idx};  
                    m_break_folding = L{5,molecule_idx}; 
                    [new_molecules]=break_molecule(c{mol,2},m_break,m_break_OCR,m_break_bl,m_break_E); % New molecules include the folding, of course.
                    M_to_be_added=[M_to_be_added new_molecules];            
    
                    HYDROLYZED(end+1,:)={m_break, m_break_OCR, m_break_bl, m_break_E, m_break_folding, new_molecules, i}; % Every row stores a molecule, and the columns mean: 
                       % "broken molecule sequence" "OCR_label" "bl_label" "E_label" "folding" {cell array with molecules formed after hydrolysis} "time step"
                end
    
                % Substract the hydrolyzed molecules from the pool:
                L(:,M_to_be_removed)=[];
    
                % Include new molecules (obtained from hydrolysis) to the pool:
                L = [M_to_be_added, L];
            end
        end 
        
        %% EVALUATE RUPTURE PROBABILITIES OF ALL THE MOLECULES IN LEVEL 0, 1 OR 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -- The order is not important, so we are going to start testing level 1 and 2 (pool) linkages (hydrogen bonds)
        
        % LEVEL 2 (POOL) => Evaluate the denaturation probability of all dsRNA AND FOLDED ssRNA molecules present in the pool in the current time step
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
                                    is_R=0;
                                    is_C=0;
                                    if Rup==1                            
                                        noligo_r=oligo_r(1:k(1)-1); % Delete possible empty ('-') positions of right oligo
                                        nOCR_r=OCR_r(1:k(1)-1);
                                        nE_r=E_r(1:k(1)-1);
                                        if any(nOCR_r=='R')
                                            R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r, []}; 
                                            is_R=1;
                                        end

                                        if any(OCR_r=='C')
                                            C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r, []};
                                            is_C=1;
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
                                        L=ampliation_E1_MODIF(L,noligo_r,new_OCR_r,bl_r,nE_r); % Introducing folding inside ampliation function
                                        % -- Copy the folding into C/R_CELL too:
                                        if is_R
                                            R_CELL{end,5}=L{5,1};
                                        end
                                        if is_C
                                            C_CELL{end,5}=L{5,1};
                                        end

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
                                        is_R=0;
                                        is_C=0;

                                        if any(nOCR_r=='R')
                                            R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r, []};  
                                            is_R=1;
                                        end

                                        if any(OCR_r=='C')
                                            C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r, []};
                                            is_C=1;
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
                                        % -- Copy the folding into C/R_CELL too:
                                        if is_R
                                            R_CELL{end,5}=L{5,1};
                                        end
                                        if is_C
                                            C_CELL{end,5}=L{5,1};
                                        end
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

                                        is_R=0;
                                        is_C=0;

                                        if any(nOCR_r=='R')
                                            R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r, []};    
                                            is_R=1;
                                        end

                                        if any(OCR_r=='C')
                                            C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r, []};
                                            is_C=1;
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
                                        % -- Copy the folding into R_CELL too:
                                        if is_R
                                            R_CELL{end,5}=L{5,1};
                                        end
                                        if is_C
                                            C_CELL{end,5}=L{5,1};
                                        end
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

                                        is_R=0;
                                        is_C=0;

                                        if any(nOCR_r=='R')
                                            R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r, []};  
                                            is_R=1;
                                        end

                                        if any(OCR_r=='C')
                                            C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l((k(z)+1):(k(z+1)-1)),nE_r, []};
                                            is_C=1;
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
                                        % -- Copy the folding into C/R_CELL too:
                                        if is_R
                                            R_CELL{end,5}=L{5,1};
                                        end
                                        if is_C
                                            C_CELL{end,5}=L{5,1};
                                        end
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

                                is_R=0;
                                is_C=0;

                                if any(nOCR_r=='R')
                                    R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r, []};     
                                    is_R=1;
                                end

                                if any(OCR_r=='C')
                                    C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(1:k(1)-1),nE_r, []};
                                    is_C=1;
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
                                % -- Copy the folding into C/R_CELL too:
                                if is_R
                                    R_CELL{end,5}=L{5,1};
                                end
                                if is_C
                                    C_CELL{end,5}=L{5,1};
                                end
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

                                is_R=0;
                                is_C=0;

                                if any(nOCR_r=='R')
                                    R_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r, []};    
                                    is_R=1;
                                end

                                if any(OCR_r=='C')
                                    C_CELL(end+1,:)={nOCR_r, noligo_r, bl_l(k(temps)+1:size(oligo_l,2)),nE_r, []};
                                    is_C=1;
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
                                % -- Copy the folding into C/R_CELL too:
                                if is_R
                                    R_CELL{end,5}=L{5,1};
                                end
                                if is_C
                                    C_CELL{end,5}=L{5,1};
                                end
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
                            is_R=0;
                            is_C=0;
                            if size(nE_r,2)==1
                                nE_r='!';
                            end

                            if any(OCR_r=='R')
                                R_CELL(end+1,:)={OCR_r, oligo_r, bl_l, E_r, []}; 
                                is_R=1;
                            end
                            
                            if any(OCR_r=='C')
                                C_CELL(end+1,:)={OCR_r, oligo_r, bl_l, E_r, []};
                                is_C=1;
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
                            % -- Copy the folding into C/R_CELL too:
                            if is_R
                                R_CELL{end,5}=L{5,1};
                            end
                            if is_C
                                C_CELL{end,5}=L{5,1};
                            end
                            
                        end
                    end
                end
            
            else

                f_mol=L{5,pos};
                f_state=f_mol{3};

                if f_state
                    f_seq=f_mol{1};
                    num_pb=sum(f_seq=='(');
                    Rup = rupture(beta(i),num_pb);
                    if Rup
                        L{5,pos}={f_mol{1},f_mol{2},[0]};
                        total_denat_fold=total_denat_fold+1;
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

                            is_R=0;
                            is_C=0;

                            if any(OCR_chain=='R')
                                R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain, []};      
                                is_R=1;
                            end

                            if any(OCR_chain=='C')
                                C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain, []};
                                is_C=1;
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
                            % -- Copy the folding into C/R_CELL too:
                            if is_R
                                R_CELL{end,5}=L{5,1};
                            end
                            if is_C
                                C_CELL{end,5}=L{5,1};
                            end
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

                            is_R=0;
                            is_C=0;

                            if any(OCR_chain=='R')
                                R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten,E_chain, []}; 
                                is_R=1;
                            end

                            if any(OCR_chain=='C')
                                C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten,E_chain, []};
                                is_C=1;
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
                            % -- Copy the folding into C/R_CELL too:
                            if is_R
                                R_CELL{end,5}=L{5,1};
                            end
                            if is_C
                                C_CELL{end,5}=L{5,1};
                            end
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

                            is_R=0;
                            is_C=0;

                            if any(OCR_chain=='R')
                                R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain, []};       
                                is_R=1;
                            end

                            if any(OCR_chain=='C')
                                C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain, []};
                                is_C=1;
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
                            % -- Copy the folding into C/R_CELL too:
                            if is_R
                                R_CELL{end,5}=L{5,1};
                            end
                            if is_C
                                C_CELL{end,5}=L{5,1};
                            end
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
                            
                            is_R=0;
                            is_C=0;

                            if any(OCR_chain=='R')
                                R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain, []};   
                                is_R=1;
                            end

                            if any(OCR_chain=='C')
                                C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain, []};
                                is_C=1;
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
                            % -- Copy the folding into C/R_CELL too:
                            if is_R
                                R_CELL{end,5}=L{5,1};
                            end
                            if is_C
                                C_CELL{end,5}=L{5,1};
                            end
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

                    is_R=0;
                    is_C=0;

                    if any(OCR_chain=='R')
                        R_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain, []};  
                        is_R=1;
                    end

                    if any(OCR_chain=='C')
                        C_CELL(end+1,:)={OCR_chain, nt_chain, bl_rewritten, E_chain, []};
                        is_C=1;
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
                    % -- Copy the folding into C/R_CELL too:
                    if is_R
                        R_CELL{end,5}=L{5,1};
                    end
                    if is_C
                        C_CELL{end,5}=L{5,1};
                    end
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

                    is_R=0;
                    is_C=0;

                    if any(OCR_chain=='R')
                        R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain, []};
                        is_R=1;
                    end

                    if any(OCR_chain=='C')
                        C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain, []};
                        is_C=1;
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
                    % -- Copy the folding into C/R_CELL too:
                    if is_R
                        R_CELL{end,5}=L{5,1};
                    end
                    if is_C
                        C_CELL{end,5}=L{5,1};
                    end
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

                is_R=0;
                is_C=0;

                if any(OCR_chain=='R')
                    R_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain, []};    
                    is_R=1;
                end

                if any(OCR_chain=='C')
                    C_CELL(end+1,:)={OCR_chain, nt_chain, bl_chain, E_chain, []};
                    is_C=1;
                end
                

                bl_new='';
                bl_new(1:size(nt_chain,2))='0'; % All links are "good", because they were generated by template depending polymerization
                bl_new(end)='X';
                new_OCR_chain=OCR_chain;
                new_OCR_chain(OCR_chain=='R')='O';   

                L=ampliation_E1_MODIF(L,nt_chain,new_OCR_chain,bl_new,E_chain);
                % -- Copy the folding into C/R_CELL too:
                if is_R
                    R_CELL{end,5}=L{5,1};
                end
                if is_C
                    C_CELL{end,5}=L{5,1};
                end
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

% ---- Update variables (more efficient code, although it increases ~30
% seconds the simulation time and the weight of the file goes from 0.2 MB to 11 MB):.
    pol_cols=cellfun(@(x) length(x)>1, L(1,:)); % <--- Just save the polymers
    seqs_data=L(1,pol_cols);
    Lt{i,1}=seqs_data; 
    % Lt_OCR{i,1}=L(2,pol_cols);
    % Lt_bl{i,1}=L(3,pol_cols);
    % Lt_E{i,1}=L(4,pol_cols);
    % Lt_fold{i,1}=L(5,pol_cols); % <-- Esta info realmente no es necesaria, teniendo la secuencia la puedo calcular
    structures=L(5,pol_cols);
    folded_fout=cellfun(@(x) x{3}==1, structures);
    folded_how_many{i,1}=sum(folded_fout); % Stores only the ssRNAs that are folded in this time step (so if they have mfe but where denaturalized in preivous step, they are not included)
    folded_how_many{i,2}=cellfun(@length,seqs_data(folded_fout)); % Stores the length of those polymers.
    total_pols_pool{i,1}=size(structures,2);
    total_pols_pool{i,2}=cellfun(@length,seqs_data);
    mfe_fout=cellfun(@(x) any(x{1}~='.'), structures);
    mfe_how_many{i,1}=sum(mfe_fout);
    mfe_how_many{i,2}=cellfun(@length,seqs_data(mfe_fout));
    mean_energy(i,1)=mean(cellfun(@(x) x{2}, structures(mfe_fout)));

   % Number and size of molecules attached to clay or level 1:
    d = diff([0 O_C 0]);
    start_idx = find(d == 1);
    end_idx   = find(d == -1) - 1;
    poly_sizes = end_idx - start_idx + 1;
    poly_sizes = poly_sizes(poly_sizes>1);
    num_polymers = numel(poly_sizes);
    total_pols_C{i,1} = num_polymers;
    total_pols_C{i,2} = poly_sizes;

    d = diff([0 O_U 0]);
    start_idx = find(d == 1);
    end_idx   = find(d == -1) - 1;
    poly_sizes = end_idx - start_idx + 1;
    poly_sizes = poly_sizes(poly_sizes>1);
    num_polymers = numel(poly_sizes);
    total_pols_U{i,1} = num_polymers;
    total_pols_U{i,2} = poly_sizes;


% ---- Update variables. NOTICE: Uncomment following lines to obtain all possible data from simulations in addition to the Activity, R_CELL and C_CELL variables (althoug notice that this will require additional time per simulation)----
 
% Ot_C=[Ot_C;O_C];
% Ot_U=[Ot_U;O_U];
% Ct=[Ct;C];
% Ct_OCR=[Ct_OCR;C_OCR];
% Ct_bl=[Ct_bl;C_bl];
% Ct_E=[Ct_E;C_E];
% Ut=[Ut;U];
% Ut_OCR=[Ut_OCR;U_OCR];
% Ut_E=[Ut_E;U_E];

         tf=i; % Utilitarian variable

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
