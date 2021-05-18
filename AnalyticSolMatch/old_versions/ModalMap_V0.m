% First code to check analytical calculations: checks for the solutions in
% matrix form. We find that the matrix-form analytical solutions match with
% the direct numerical solutions

% @pragya: Fall 2019

close all;
clear all; 
% import python module
py.importlib.import_module('networkx');

% ------------------ SPECIFY PARAMETERS ----------------------------------
% Network parameters -----------------------------------------------------
% number of nodes
N = 50; 
% probability of rewiring (for WS graph)
p = 0.1; 
% set intra-layer directedness
directed = 0;  
% set the magnitude of interlayer connections
dK = 1; 
% interlayer network 
kappa = dK*eye(N); 

% Dynamic parameters ------------------------------------------------------ 
% total time
T = 1; 
nt = 100; 
tarray =  linspace(0,T,nt+1); 
delt = T/nt ; 
Niter = 100; 

% -------------------------------------------------------------------------
% ---------------------- CONSTRUCT DUPLEX ---------------------------------
% Create duplex : layer 1 is WS and layer 2 is BA of comparable density ---

G_WS = net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(...
                                                 int16(N),int16(N/2),0.1));
                                             
G_BA = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(...
                                                     int16(N),int16(N/2)));
%G_BA = G_WS;                                                  
Duplex0 = zeros(2*N, 2*N); 
Duplex0(1:N,1:N)= G_WS; Duplex0(N+1: 2*N, N+1: 2*N) = G_WS; 

% symmetrize only intra-layer networks
Duplex0 = Duplex0 + Duplex0' - Duplex0.*eye(size(Duplex0,1)); 

% set off-diagonal matrix
Duplex0(N+1:2*N,1:N) = kappa; 

% stabilize the network
Duplex = Duplex0/(eigs(Duplex0,1)+1) - eye(size(Duplex0));
[VL1, DL1] =  eig(Duplex(1:N,1:N));  vec = diag(DL1); 
norm_max = max(max(abs(vec))); 

% normalize by most negative (largest in magnitude) eigenvalue of layer 1
Duplex = Duplex/norm_max; 

% extract stabilized layer networks
Sym_WS = Duplex(1:N,1:N); 
Sym_BA = Duplex(N+1:2*N, N+1:2*N); 
off_mat = Duplex(N+1:end,1:N); 

% extract scaling of off-diagonal connections  
off_diag = diag(off_mat); 
dscale = off_diag(1);  

% ------------------- Get network properties ------------------------------
[V1,D1] = eig(Sym_BA);      [V2,D2] = eig(Sym_WS); 

% sort the eigenvalues in order of decreasing absolute value of the eigenvalues
% such that most negative eigenvalue (largest magnitude) is the first
[D1,I1] = sort(diag(D1));  V1 = V1(:,I1); 
[D2,I2] = sort(diag(D2));  V2 = V2(:,I2); 
muvec = D2;           xivec = D1; 

% calculate spectral distance and alignment matrix; 
dspec = sqrt(sum((D1- D2).^2)); 
C_mat = zeros(N,N); 
for i =1:N
    for j=1:N 
      C_mat(i,j) = V2(:,i)'*V1(:,j); 
    end
end



XF = [0,1]; 

for ixi =1:N
    for imu =1:N
        Clm = C_mat(ixi,imu); 
        [state,costate,E(ixi,imu)] =...
                        OneModeSol_funV0(xivec(ixi), muvec(imu),Clm, T, nt, XF); 
        
    end
end

imagesc(E)








