function [H_star, F, A_star,S_opt,W_opt,C_opt, Obj] = IGEEGc(A_x, A_h,M, parameter)
num_sample = parameter.n;
num_cluster = parameter.c;
tol = 1e-6;
max_iter = 30;
tol_iter = 2;
lambda = parameter.lambda;
%% initialize H_star F C S
A_star = (1-lambda)*A_x(:,:)+lambda*A_h(:,:);
[H_star,C_0,S_0, F_0, L] = initialHFCS(A_star,M, parameter);
F = F_0;
S_tp = S_0;
W_tp = parameter.W;
C_tp = C_0;
[L_S_tp,~]=lapMatrix(S_tp);
%% loop
Obj = zeros(1,30);
for iter = 1 : max_iter
     Obj(iter) = obj_function(H_star,A_star,F,S_tp,W_tp,C_tp,L_S_tp,M,parameter);
    if iter > tol_iter &&  abs((Obj(iter) - Obj(iter-1) )/Obj(iter-1))< tol
        break;
    end    
    % update S* W C
    [L_S_tp S_tp] = update_S(S_0,M,W_tp,C_tp,parameter);
    W_tp = update_W(S_tp,M,parameter);
    C_tp = update_C(C_0,S_tp,parameter);
    % update A
    % A*=（1-q）x Ax + q x Ah 
    [H, ~, A_x] = spectral_embedding(S_tp, num_cluster);
    X = sqrt(sum(H.^2, 2));
    Hn = H ./ repmat(X, 1, size(H, 2));
    affn = constructW_Gaussian(Hn',parameter.sigma);
    [~, ~, A_h] = spectral_embedding(affn, num_cluster);
    A_star = (1-lambda)*A_x(:,:)+lambda*A_h(:,:);
    % update H
     [Uh, ~, Vh] = svd(A_star*F, 'econ');
     H_star = Uh * Vh';
   
    % update F
     G = A_star'*H_star;
    L = update_L(L,G, num_sample, num_cluster);
    F = L*(L'*L + eps*eye(num_cluster))^(-0.5);
    
end
S_opt = S_tp;
W_opt = W_tp;
C_opt = C_tp;
disp(['inter iter:', num2str(iter-1)]);
end
