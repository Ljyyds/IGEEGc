function Obj=obj_function(H,A,F,S,W,C,L_S,M,parameter)
So = S-twodi23D(M,W);
part1=0.5*parameter.gamma1*trace(So'*So);
part2=0.5*parameter.gamma2*trace(C'*L_S*C);
part3=trace(H' * A * F);
Obj=part1+part2-part3;
end