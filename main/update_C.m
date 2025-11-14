function C_tp1=update_C(C_0,S_t,parameter)

[mS nS]=size(S_t);
[Lap D]=lapMatrix(S_t);
C_t=C_0;

t=1;
while t<parameter.iter+1
    derR=parameter.gamma2*Lap*C_t; % derivative of R
    C_t=C_t-parameter.eta*derR;
    t=t+1;
end
C_tp1=C_t;
end





