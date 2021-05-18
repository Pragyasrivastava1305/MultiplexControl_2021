% This code matches solutions for costate variables.
% Initial conditions are obtained by calling the function min_eng_cont.
% Variables are looped over rather than being in vectoised form. 
% Given correct initial conditions, this code checks for a correct
%  dynamical evolution of costate variables.
% pragyasr : Fall 2019

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
                                                 
Duplex0 = zeros(2*N, 2*N); 
Duplex0(1:N,1:N)= G_WS; Duplex0(N+1: 2*N, N+1: 2*N) = G_BA; 

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
[V1,D1] = eig(Sym_WS);      [V2,D2] = eig(Sym_BA); 

% sort the eigenvalues in order of decreasing absolute value of the eigenvalues
% such that most negative eigenvalue (largest magnitude) is the first
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


% ------------------- ALLOCATIONS/INITIALIZATION --------------------------
% state vectors 
x1sol = zeros(N,nt+1);  x2sol = zeros(N,nt+1); 

% costate vectors
lamsol= zeros(N,nt+1); sigsol= zeros(N,nt+1); 

% initial condition for state variables
X1_init = zeros(N,1);    X2_init = zeros(N,1);     X0 = [X1_init; X2_init];

% input matrix
B = zeros(2*N,N);       B(1:N,1:N) = eye(N); 

% -------------------- ITERATIONS OVER SOLUTIONS  -------------------------
%for iter =1: Niter
    iter =1; 
    rng(iter);     X1_final = randn(N,1);
    rng(2*iter);   X2_final = randn(N,1); 
    XF = [X1_final; X2_final];              XF = XF/norm(XF);
    
    % numerical solution --------------------------------------------------
	[x,u,v,n_err] = min_eng_cont(Duplex, T, B, X0, XF, nt,0); 
    costate_num = v(2*N+1:end,:); 
    
    % extract initial state of costate variables 
    sig_num = costate_num(N+1:end,:);
    lam_num = costate_num(1:N,:); 
    
    % calculate energy
    u = u'; 
    for ilen = 1:size(u,2)    
		E(ilen) = delt*trapz(u(:,ilen).*u(:,ilen)); 
    end
    minE_num(iter) = sum(E);

%end

% Analytical solutions  ---------------------------------------------------
L0_num = v(2*N+1:end,1);
lam0 = pinv(V1)*L0_num(1:N); 
sig0 = pinv(V2)*L0_num(N+1:end); 
muvec = D2;           xivec = D1; 

lam_inter = zeros(N,size(tarray,2)); 
sig_inter = zeros(N,size(tarray,2)); 
for iN =1:N
	sig_inter(iN,:) = sig0(iN)*exp(-muvec(iN)*tarray);

    fun_ac = zeros(1,size(tarray,2)); 
    for jN = 1:N
        fun_ac  = fun_ac +  dscale*C_mat(jN,iN) * sig0(jN)...
                                                * (exp(-muvec(jN)*tarray)...
                                                   - exp(-xivec(iN)*tarray))...
                                                /(muvec(jN) - xivec(iN)); 
    end
    lam_inter(iN,:) = exp(-xivec(iN)*tarray)*lam0(iN) +  fun_ac;
end

% transform back to original basis ----------------------------------------
lam_comp = zeros(N,nt+1); 
sig_comp = zeros(N,nt+1); 

for it = 1:nt+1
    sig_comp(:,it) = V2*sig_inter(:,it); 
    lam_comp(:,it) = V1*lam_inter(:,it); 
end
sup_mat = [lam_comp; sig_comp]; 

% make plots to compare ---------------------------------------------------
imagesc(round(sup_mat - v(2*N+1:end,:),6)); colorbar
max(max(abs(sup_mat - v(2*N+1:end,:))))




































