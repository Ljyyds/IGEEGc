function X=twodi23D(X_ori,C)
%三维转二维
%C
%X_ori
[t,n,m]=size(X_ori);
X=zeros(n,m);
for i=1:t
    X_slice=C(i).*X_ori(i,:,:);
    X_slice=reshape(X_slice,n,m);
    X=X+X_slice;
end
end