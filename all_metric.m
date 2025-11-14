function [Result] = all_metric(L_opt,L_true,parameter)

NMI_1=[];
NMI_2=[];
SR_t=[];
SAR_t=[];
microfscore=[];
ka = [];
for i=1:parameter.c
    for j=1:parameter.c
        NMI_1=[NMI_1 NMI(L_opt(i,:),L_true(j,:))];
        [zRand,SR,SAR,VI]=zrand(L_opt(i,:),L_true(j,:));
        SR_t=[SR_t SR];
        SAR_t=[SAR_t SAR];
        [micro, macro]=micro_macro_PR(L_opt(i,:),L_true(j,:));
        microfscore=[microfscore micro.fscore];


        kappa= kappaindex(L_opt(1,:)+1,L_true(1,:)+1,parameter.c);
        %---------------
        kappa= kappaindex(L_opt(i,:)+1,L_true(j,:)+1,parameter.c);
        ka=[ka kappa];
        microfscore=max(microfscore);
    end
end
NMI_1=max(NMI_1);
SAR_F=max(SAR_t);
microfscore=max(microfscore);
kappa= max(ka);
RI = RandIndex(L_opt,L_true); % Calculate Rand Index
Result=[RI NMI_1 SAR_F microfscore kappa ........
    parameter.sigma parameter.lambda parameter.gamma1 parameter.gamma2 parameter.gamma3];
end