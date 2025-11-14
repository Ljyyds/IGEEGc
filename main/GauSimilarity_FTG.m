function [M M_2D W]=GauSimilarity_FTG(FTGs,parameter)
    W = parameter.W;
    sigma = parameter.sigma;
    mE = length(FTGs);
    D_S = zeros(3,mE,mE);

    dists = zeros(3,mE,mE);

    for t = 1:3
        disti = zeros(mE,mE);
        parfor i = 1:mE
            FTGi = FTGs{i};
            temp = zeros(1,mE);
            for j = i:mE
                if i == j 
                    temp(j) = 0;
                    continue;
                end
                FTGj = FTGs{j};
                temp(j) = DTW_FTG(FTGi,FTGj,t);
            end
            disti(i,:) = temp;
        end
        dists(t,:,:) = triu(disti) + triu(disti,1).';
    end

distg = z_regularization(reshape(dists(1,:,:),mE,mE));
distk = z_regularization(reshape(dists(2,:,:),mE,mE));
distw = z_regularization(reshape(dists(3,:,:),mE,mE));

dists_z(1,:,:)=distg;
dists_z(3,:,:)=distk;
dists_z(2,:,:)=distw;
M = zeros(3,mE,mE);

for t = 1:3
    for i = 1:mE
        A = reshape(dists_z(t,i,:),[1,mE]);
        for j = i:mE
            B = reshape(dists_z(t,j,:),[1,mE]);
            g = norm(A-B);
            M(t,i,j) = exp(-g^2/(2*sigma^2));
            M(t,j,i) = M(t,i,j);
        end
    end
end

M_2D  = twodi23D(M,W);
end

