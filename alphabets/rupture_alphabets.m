function [Rup]=rupture_alphabets(beta,R1,R2,alphabet)
% Beta Emviormental conditions for rupture
% R1 is the strand we want to address the number of hydrogen bonds and R2
% is its complementary strand.
Sc=size(R1,2);
% Sc is the size of the chain

if strcmp(alphabet,'2_2') % alphabet 2, with A-U baise pair (2 hydrogen bonds)
    k=2*Sc;
elseif strcmp(alphabet,'2_3') % alphabet 2, with G-C baise pair (3 hydrogen bonds)
    k=3*Sc;
elseif strcmp(alphabet,'4') % common alphabet 4 => A-U (2pdH), G-C (3pdH)
    num_g=count(R1,'g');
    num_c=count(R1,'c');
    k=num_g+num_c+(2*Sc);
elseif strcmp(alphabet,'4*') % four letters, 3 interationcs => A-U (2pdH), G-C (3pdH), G-U (2pdH)
    num_c1=count(R1,'c');
    num_c2=count(R2,'c');
    k=num_c1+num_c2+(2*Sc);
elseif strcmp(alphabet,'6_2') % alphabet 6, with xy baise pair displaying 2 hydrogen bonds
    num_g=count(R1,'g');
    num_c=count(R1,'c');
    k=num_g+num_c+(2*Sc);
elseif strcmp(alphabet,'6_3') % alphabet 6, with xy baise pair displaying 3 hydrogen bonds
    num_a=count(R1,'a');
    num_u=count(R1,'u');
    k=(3*Sc)-num_a-num_u;
else
    disp('No such alphabet, please check your initial conditions!')
end

ProbRup = exp(-beta*k); % Probability of rupture 

na=rand;

if na < ProbRup
    Rup=1;
else
    Rup=0;
end