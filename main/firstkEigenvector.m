function [firstkVector,firstkLambda]=firstkEigenvector(Ei,k) 
%First k eigenvectors corresponding to k smallest eigenvalues

matrixSize=size(Ei);
[v,d] = eig(Ei);
v = [v ; diag(d)'; abs(diag(d))']';
v = sortrows(v,matrixSize+2); %Eigenvalues in ascending order
firstkLambda = v(1:k,matrixSize(1)+1) ;
firstkVector = v(1:k,1:matrixSize)';