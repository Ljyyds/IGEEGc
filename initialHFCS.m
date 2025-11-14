function [H_star,C_0,S_0, F, L] = initialHFCS(A_star,M, parameter)
%% initialize H_star
num_cluster = parameter.c;
[H_star] = eig1(A_star, num_cluster);
%% initialize L F
stream = RandStream.getGlobalStream;
reset(stream);
X = sqrt(sum(H_star .^ 2, 2));
H_normalized = H_star ./ repmat(X, 1, size(H_star, 2));
label = kmeans(H_normalized, num_cluster, 'maxiter', 1000, 'replicates', 20, 'emptyaction', 'singleton');
L = idx2pm(label);
F = L * (L' * L + eps * eye(num_cluster)) ^ ( - 0.5);
%% initialize S0
mD = parameter.n;
S_0=1/mD*ones(mD);
%% initialize C0
sen=zeros(mD,1);   
for t=1:parameter.c
    for j=1:parameter.v
        sen((mD/parameter.c*(t-1)+mD/(parameter.c*parameter.v)*(j-1)+1):(mD/parameter.c*(t-1)+mD/(parameter.c*parameter.v)*j))=j;
    end
end
M_2D = twodi23D(M,parameter.W);
degree_vector = sum(M_2D, 2);   
D_M = diag(degree_vector);
L_M = D_M-M_2D;
sen_unique=unique(sen);
v=length(sen_unique);
sen_unique=reshape(sen_unique,[1, v]);
senUpd=sen;
temp=1;
for t=sen_unique
    senUpd(sen==t)=temp;
    temp=temp+1;
end
B=zeros(mD,v);
for t=1:v-1
    temp=(senUpd==t);
    B(temp,t)=1;
    groupsize=sum(temp);
    B(:,t)=B(:,t)-groupsize/mD;
end

Q=null(B');
Z=sqrtm(Q'*D_M*Q);


invZ=inv(Z);

forU=invZ'*Q'*L_M*Q*invZ;
forU=(forU+forU')/2;

[U,firstkLambda]=firstkEigenvector(forU,parameter.c);

C_0=Q*invZ*U;
end



