function [A_x, A_h,M,W] = PreProcess(FTGs, parameter)
num_cluster = parameter.c;
[M,M_2D,W] = GauSimilarity_FTG(FTGs,parameter); 
[H, ~, A_x] = spectral_embedding(M_2D, num_cluster);
A_x = full(A_x);
X = sqrt(sum(H.^2, 2));
Hn = H ./ repmat(X, 1, size(H, 2));
affn = constructW_Gaussian(Hn',parameter.sigma);
[~, ~, A_h] = spectral_embedding(affn, num_cluster);
end

