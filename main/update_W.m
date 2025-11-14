function [W_i]=update_W(S_tp1,M_0,parameter)

[~, n, ~] = size(M_0);
k = 1;
W = zeros(n * (n + 1) / 2, size(M_0, 1));
random = parameter.nn;

index_map = zeros(n * (n + 1) / 2, 2);
k = 1;

for i = 1:n
    for j = i:n
        index_map(k, :) = [i, j];
        k = k + 1;
    end
end

parfor idx = 1:size(index_map, 1)
    i = index_map(idx, 1);
    j = index_map(idx, 2);
    S_col = M_0(:, random(i), random(j));
    pinvS = pinv(S_col);
    W(idx, :) = S_tp1(random(i), random(j)) .* pinvS; % 逐元素乘法
end

[nC,~]=size(W);
W_i=sum(W,1)/nC;
W_i=W_i/sum(W_i);

end

