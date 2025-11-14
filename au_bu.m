function [a_opt,b_opt] = au_bu(U,i)
%A 此处显示有关此函数的摘要
%   此处显示详细说明

rep = mean(U{i});
%rep = U{i}(length(U{i})/2);
amin = min(U{i});
amax = rep;
bmin = rep;
bmax = max(U{i});

fa = @(au)-(sum(U{i}>au&U{i}<rep)/sum(U{i}<rep))*(1-(rep-au)/(rep-min(U{i})));
fb = @(bu)-(sum(U{i}>rep&U{i}<bu)/sum(U{i}>rep))*(1-(bu-rep)/(max(U{i})-rep));

fau = @(au)-(sum(U{i}>au&U{i}<rep))*exp(-(rep-au));
fbu = @(bu)-(sum(U{i}<bu&U{i}>rep))*exp(-(bu-rep));

aa = 1;
fa1 = @(au)-((sum(U{i}>au&U{i}<rep)/sum(U{i}<rep)))*exp(-((rep-au)/(max(U{i})-min(U{i}))).^aa);
fb1 = @(bu)-((sum(U{i}<bu&U{i}>rep)/sum(U{i}>rep)))*exp(-((bu-rep)/(max(U{i})-min(U{i}))).^aa);

fa2 = @(au)-sum(U{i}>au&U{i}<rep)*exp(-(aa*(rep-au)));
fb2 = @(bu)-sum(U{i}<bu&U{i}>rep)*exp(-(aa*(bu-rep)));
[a_opt,~] = fminbnd(fa1,amin,amax);
[b_opt,~] = fminbnd(fb1,bmin,bmax);

end

