close all;
clear all;
clc;


name = ['Data_IV_3_s2'];
load([name '.mat']);
data=eval(name);
[mD nD]=size(data);
parameter.c=max(data(:,1)); 
parameter.n = mD;            
parameter.fs = 400;  %173.61 bonn         
parameter.ftigs = 20;  %object 
parameter.y = 0.01;  %l1 trend filter parameter

FTGs = granulations(data(:,2:end),parameter);
num = parameter.ftigs;
save(name,'FTGs',name,'num');
