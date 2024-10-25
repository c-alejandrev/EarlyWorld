function [Rup]=rupture(beta,Sc)
% j= Beta Emviormental conditions for rupture
% Sc is the size of the chain

ProbRup = exp(-beta*Sc); % Probability of rupture 
% rng shuffle
% Decide if the chain will break or not
% Rup = binornd(1,ProbRup); 
% Rup=sum(rand >= cumsum([0,1-ProbRup,ProbRup]))-1;

na=rand;

if na < ProbRup
    Rup=1;
else
    Rup=0;
end

% if Sc > 5 
% %     if Rup==0
% %     0
% %     end
%     if Rup==1
% Sc
% ProbRup
% na
%     end
% end

