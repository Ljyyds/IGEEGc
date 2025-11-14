function tg_dist = DTW_FTG(FTG1, FTG2, num)
    M = length(FTG1.FTIGs);
    N = length(FTG2.FTIGs);

    a1 = [FTG1.FTIGs.aopt]; b1 = [FTG1.FTIGs.bopt];
    ki1 = [FTG1.FTIGs.ki]; W1 = [FTG1.FTIGs.W];
    a2 = [FTG2.FTIGs.aopt]; b2 = [FTG2.FTIGs.bopt];
    ki2 = [FTG2.FTIGs.ki]; W2 = [FTG2.FTIGs.W];
    bps1 = FTG1.bps; bps2 = FTG2.bps;

    if isequal(FTG1, FTG2)
        tg_dist = 0;
        return;
    end

    costMatrix = zeros(M, N);

    costMatrix(1,1) = dist(a1,b1,ki1,W1,bps1,a2,b2,ki2,W2,bps2,1,1,num);
    for i = 2:M
        costMatrix(i,1) = costMatrix(i-1,1) + ...
            dist(a1,b1,ki1,W1,bps1,a2,b2,ki2,W2,bps2,i,1,num);
    end
    for j = 2:N
        costMatrix(1,j) = costMatrix(1,j-1) + ...
            dist(a1,b1,ki1,W1,bps1,a2,b2,ki2,W2,bps2,1,j,num);
    end

    for i = 2:M
        for j = 2:N
            d = dist(a1,b1,ki1,W1,bps1,a2,b2,ki2,W2,bps2,i,j,num);
            costMatrix(i,j) = min([costMatrix(i-1,j-1), costMatrix(i,j-1), costMatrix(i-1,j)]) + d;
        end
    end

    tg_dist = costMatrix(M, N);
end

function d = dist(a1,b1,ki1,W1,bps1,a2,b2,ki2,W2,bps2,i,j,num)
    switch num
        case 1
            d = max(abs(a1(i)-a2(j)), abs(b1(i)-b2(j)));
        case 2
            d = ((bps1(i+1)-bps1(i)) + (bps2(j+1)-bps2(j))) * abs(ki1(i)-ki2(j)) / 2;
        case 3
            d = abs(W1(i)-W2(j));
    end
end