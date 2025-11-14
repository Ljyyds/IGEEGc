clear all;
close all;
clc;
warning off;

path    = 'dataset/';
dataset = 'Data_IV_2a_s1';
load ([path dataset,'.mat']);

data   = eval(dataset);
L_true = data(:,1); % true labels without labled EEG   
data   = data(:,2:end);
parameter.c = max(L_true); % clusting num
parameter.v = parameter.c; % fairness (label ratio)
parameter.n = size(data,1); % num of data
parameter.W = [1/3,1/3,1/3];

parameter.sigma  = 1;
parameter.lambda = 1;
parameter.gamma1 = 1;
parameter.gamma2 = 1;
parameter.gamma3 = 1;

parameter.iter = 50;
parameter.eta  = 0.001;
parameter.nn   = randperm(parameter.n);% Random disruption is used to update W


[A_x, A_h,M,parameter.W] = PreProcess(FTGs, parameter);% compute A_x A_h and the similarity matrix M
[H_star, F, A_star,S_opt,W_opt,C_opt, Obj] = IGEEGc(A_x, A_h,M, parameter);


L_opt = Lfinal(F);
L_true = l_true_reshape(L_true,parameter.c);

Result = all_metric(L_opt,L_true,parameter);

