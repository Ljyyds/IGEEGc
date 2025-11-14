function [F] = Lfinal(F)
F = F';
[mL,nL]=size(F); 
for i=1:nL
    L_max=max(F(:,i));
    for j=1:mL
        if F(j,i)==L_max
            F(j,i)=1;
        else
            F(j,i)=0;
        end
    end
end
end