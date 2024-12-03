function [Rup]=rupture_alphabets(beta,R1,R2,alphabet)
% beta stands for the value of beta
% R1 is the strand we want to address the number of hydrogen bonds and R2 is its complementary strand.
% alphabet is a string code that represents a particular alphabet (here six different alphabets are studied, but more can be added)

Sc=size(R1,2);
% Sc is the length of the molecule

% In the following if-else loop we calculate the total number of hydrogen bonds (k) that will depend on the alphabet we use:
if strcmp(alphabet,'2_2') % alphabet of 2 letters, with A-U baise pair (2 hydrogen bonds)
    k=2*Sc; 
elseif strcmp(alphabet,'2_3') % alphabet of 2 letters, with G-C nucleobase pair (3 hydrogen bonds)
    k=3*Sc;
elseif strcmp(alphabet,'4') % alphabet of 4 letters => A-U (2 hydrogen bonds), G-C (3 hydrogen bonds)
    num_g=count(R1,'g');
    num_c=count(R1,'c');
    k=num_g+num_c+(2*Sc);
elseif strcmp(alphabet,'4*') % alphabet of 4 letters => A-U (2 hydrogen bonds), G-C (3 hydrogen bonds), G-U (2 hydrogen bonds)
    num_c1=count(R1,'c');
    num_c2=count(R2,'c');
    k=num_c1+num_c2+(2*Sc);
elseif strcmp(alphabet,'6_2') % alphabet of 6 letters =>  A-U (2 hydrogen bonds), G-C (3 hydrogen bonds), X-Y (2 hydrogen bonds)
    num_g=count(R1,'g');
    num_c=count(R1,'c');
    k=num_g+num_c+(2*Sc);
elseif strcmp(alphabet,'6_3') % alphabet of 6 letters => A-U (2 hydrogen bonds), G-C (3 hydrogen bonds), X-Y (3 hydrogen bonds)
    num_a=count(R1,'a');
    num_u=count(R1,'u');
    k=(3*Sc)-num_a-num_u;
else
    disp('No such alphabet, please check your initial conditions!')
end

ProbRup = exp(-beta*k); % Denaturation probability dependent on the number of hydrogen bonds (k) of the duplex R1-R2.

na=rand;

if na < ProbRup
    Rup=1;
else
    Rup=0;
end
