function [L_S_tp S_tp1] = update_S(S_0,M,W_t,C,parameter)

S_t=S_0;
[mS nS]=size(S_t);
t=1;
tmp = twodi23D(M,W_t);
P1 = S_0 - tmp;
P1 = parameter.gamma1*P1;
P2 = 0.5*parameter.gamma2*(C*C');
derS=parameter.eta*(P1-P2); % derivative of S
while t<parameter.iter+1
    S_t=S_t-derS;
    t=t+1;
end

for i=1:mS
    for j=1:mS
        if S_t(i,j)<=0
            S_t(i,j)=0;
        end
    end
end
[L_S_tp,~]=lapMatrix(S_t);
S_tp1=S_t;
end