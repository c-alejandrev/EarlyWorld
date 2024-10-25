function [Lnew]=reduction_E1(L,P1)
% If one nt or one oligomer dissapears from de solution, the nts that were
% on its right colums move one column to the left, so now we have a new L
%(Lnew) cell (matrix) of 1 unit less of column size than the original L.
s=size(L,2);

if P1==1
    Lnew=L(:,2:end);
elseif P1==s
    Lnew=L(:,1:end-1);
else
    Lnew=L(:,[1:P1-1,P1+1:end]);
end
