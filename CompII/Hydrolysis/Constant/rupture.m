function [Rup]=rupture(beta,Sc)
% beta stands for the value of alpha or beta
% Sc is the length of the molecule

ProbRup = exp(-beta*Sc); % Probability of desorption/denaturation

na=rand;

if na < ProbRup
    Rup=1;
else
    Rup=0;
end


