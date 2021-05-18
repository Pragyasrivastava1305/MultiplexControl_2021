% create stable duplex
% @pragya: spring 2020

% This code creates duplex networks. Each layer is generated from networkx
% algorithms to generate random networks. They are combined in a cuplex
% with interlayer connections given by identity matrix. 
% normalization used:  Stabilization by the real part of largest
% eigen-value of initial duplex. Then normalization by the magnitude of
% largest eigenvalue of first layer in order to measure all rates in its
% units. 

% ------------------ SPECIFY PARAMETERS ----------------------------------
% Network parameters -----------------------------------------------------
% number of nodes
%N = 3; 
% probability of rewiring (for WS graph)
p = 0.1; 
% set intra-layer directedness
directed = 0;  
% set the magnitude of interlayer connections
dK = 1; 
% interlayer network 
kappa = dK*eye(N); 
% initialize duplex
Duplex0 = zeros(2*N, 2*N); 


% Construct duplex -------------------------------------------------------
 G_WS = net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(...
                                                  int16(N),int16(2),0.1));
%                                            
   G_BA = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(...
                                                        int16(N),int16(N/2)));
%  G_BA = G_WS;
% G_WS = G_BA;                                                   
% G_latt = ones(N,N);
% G_BA = G_latt; 
% G_WS = G_BA;

Duplex0(1:N,1:N)= 0.5*(G_WS + G_WS') - G_WS.*eye(N); 
Duplex0(N+1:2*N, N+1:2*N) = 0.5*(G_BA + G_BA') - G_BA.*eye(N);  
% set off-diagonal matrix to get initial Duplex 
Duplex0(N+1:2*N,1:N) = kappa; 

% ------------------------------------------------------------------------

% stabilize the duplex
max_eig_val = real(eigs(Duplex0,1)); 
Duplex = Duplex0/(max_eig_val+1) - eye(size(Duplex0));
%Duplex = Duplex0; 
%[VL1, DL1] =  eig(Duplex(1:N,1:N));  vec = diag(DL1); 
%norm_max = max(max(abs(vec))); 
% normalize by most negative (largest in magnitude) eigenvalue of layer 1
% thus all strengths measured in terms of largest eigenvalue of layer 1
%Duplex = Duplex/norm_max; 

% extract stabilized layer networks
layer1 = Duplex(1:N,1:N); 
layer2 = Duplex(N+1:2*N, N+1:2*N); 
off_mat = Duplex(N+1:end,1:N); 

% extract scaling of off-diagonal connections  
off_diag = diag(off_mat); 
dscale = off_diag(1);  


%------------------- eigen-decomposition of each layer --------------------
[V1,D1] = eig(layer1);      
[V2,D2] = eig(layer2); 
%%%
% sort the eigenvalues in order of decreasing absolute value of the eigenvalues
% i.e. most negative eigenvalue (largest magnitude) is the first. 
% note that all eigenvalues are negative, hence default ordering of 'sort'
% function in matlab works. 

[D1,I1] = sort(diag(D1));  V1 = V1(:,I1); 
[D2,I2] = sort(diag(D2));  V2 = V2(:,I2); 

% calculate spectral distance and alignment matrix; 
dspec = sqrt(sum((D1- D2).^2)); 
C_mat = zeros(N,N); 
for i =1:N
    for j=1:N 
      C_mat(i,j) = V2(:,i)'*V1(:,j); 
    end
end






