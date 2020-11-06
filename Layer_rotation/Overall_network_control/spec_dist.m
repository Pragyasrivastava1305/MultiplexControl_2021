function [dspec,lamA,lamB] = spec_dist(A,B,opt) 
% function to calculate spectral distance
% INPUTS
% A,B   : adjacency matrices of the two networks of same size
% opt   : if opt =0, only the eigenvalues of adjacency matrix is used.
%       : if opt =1,  eigenvalues of Laplacian and adjacency matrices are
%       used
%       : if opt =2, in addition normalized Laplacian are also used 
% OUTPUTS
% dspec:  spectral distance
% lamA :  a 3N(opt==2), 2N(opt==1) or N(opt==0) length array
%         First, second and last N entries are eigenvalues 
%         of adjacency matrix, Laplacian matrix and normalized Laplacian
%         matrix for A
% lamB :  a 3N(opt==1) or 2N(opt==0)  length array. Same as lamA but for B 

N = size(A,1);

if opt ==0
    lamA = zeros(N,1); lamB = zeros(N,1); 
elseif opt == 1
    lamA = zeros(2*N,1); lamB = zeros(2*N,1); 
else
    lamA = zeros(3*N,1); lamB = zeros(3*N,1);
end
    
degA = sum(A); degB = sum(B);

DA = diag(degA); DB = diag(degB); 

LapA = DA-A; LapB = DB-B; 

lamA(1:N) = sort(eig(A),'descend');  

lamB(1:N) = sort(eig(B),'descend'); 

if opt == 1
    
    lamA(N+1:2*N) = sort(eig(LapA),'descend'); 
    
    lamB(N+1:2*N) = sort(eig(LapB),'descend'); 

elseif opt == 2
    
    sqrtDA = mpower(DA,-1/2); 
    
    sqrtDB = mpower(DB,-1/2); 
    
    normLA = sqrtDA*LapA*sqrtDA; 
    
    normLB = sqrtDB*LapB*sqrtDB;
    
    lamA(2*N+1:3*N) = sort(eig(normLA),'descend');
    
    lamB(2*N+1:3*N) = sort(eig(normLB),'descend'); 

end

dspec = sqrt(sum((lamA-lamB).^2)/length(lamA)); 



end

