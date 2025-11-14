function L = update_L(L,G, num_sample, num_cluster)
yg = diag(L'* G)';
yy = diag(L'*L+eps*eye(num_cluster))';
for i = 1 : num_sample
    gi = G(i,:);
    yi = L(i,:);
    si = (yg+gi.*(1-yi))./sqrt(yy+1-yi) - (yg-gi.*yi)./sqrt(yy-yi);
    [~,index] = max(si(:));
    L(i,:) = 0;
    L(i,index) = 1;
end
end

