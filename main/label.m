function [labels] = label(F)
[m n] = size(F);
labels = zeros(1,n);
for i = 1:n
    for j = 1:m
        if F(j,i) ==1 
            labels(i) = j;
        end
    end

end

end

