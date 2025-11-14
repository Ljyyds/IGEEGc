function [a_opt,b_opt] = au_bu(U,i)


rep = mean(U{i});
amin = min(U{i});
amax = rep;
bmin = rep;
bmax = max(U{i});

aa = 1;
fa1 = @(au)-((sum(U{i}>au&U{i}<rep)/sum(U{i}<rep)))*exp(-((rep-au)/(max(U{i})-min(U{i}))).^aa);
fb1 = @(bu)-((sum(U{i}<bu&U{i}>rep)/sum(U{i}>rep)))*exp(-((bu-rep)/(max(U{i})-min(U{i}))).^aa);


[a_opt,~] = fminbnd(fa1,amin,amax);
[b_opt,~] = fminbnd(fb1,bmin,bmax);

end

