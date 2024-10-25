function [Compl]=complementary2(D,RC1,RC2)
% Length of RC1 must be the same that RC2

RC=zeros(1,size(RC1,2));
for i=1:size(RC1,2)
    Rtemp=strcat(RC1(i),RC2(i)); % Create the pair to compare complementary with posible pairs in D
    for j=1:size(D,1)
        if Rtemp==D(j,:)
            RC(i)=1;
        elseif  fliplr(Rtemp)==D(j,:)
            RC(i)=1;
        end
    end
end

if sum(RC)==size(RC1,2)
    Compl=1;
else
    Compl=0;
end
