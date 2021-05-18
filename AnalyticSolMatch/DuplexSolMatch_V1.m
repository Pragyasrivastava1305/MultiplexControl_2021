% UPGRADES OVER V0: 
% Checks both state and costate solutions 
% Vectorized calculations instead of loop-overs
% Initial condition for costate variables still from min_eng_cont 
% Finds match: affirmative 
% @pragya: Fall 2019
% Version: 1

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
% set intra-layer directedness (layer networks are symmetric)
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
Duplex = Duplex0/(real(eigs(Duplex0,1))+1) - eye(size(Duplex0));
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



% MATRIX CHECKS -----------------------------------------------------------
Tlam = zeros(N,N,nt+1);
G = zeros(N,N,nt+1); 
H = zeros(N,N,nt+1); 
Gprime = zeros(N,N,nt+1); 
Hprime = zeros(N,N,nt+1); 
C_mat = dscale*C_mat; 

for tt = 1:nt+1
    time = tarray(tt);
    G(:,:,tt) = diag( sinh(xivec*time)./xivec/2 ); 
     
	for ii =1:N
    	for jj =1:N
       		 Tlam(ii,jj,tt) = C_mat(ii,jj)*(exp(-muvec(ii)*time) - exp(-xivec(jj)*time))...
                				                     /(muvec(ii) - xivec(jj));
                                
       		 H(ii,jj,tt) = C_mat(jj,ii)*( (exp(-muvec(jj)*time) - exp(xivec(ii)*time))...
                		                              /(muvec(jj) + xivec(ii)) + ... 
                                		        sinh(xivec(ii)*time)/xivec(ii) )/...
                                 				 (muvec(jj) - xivec(ii))/2 ;
                              
     	 	 Gprime(ii,jj,tt) = ...
       	 	  C_mat(ii,jj)*((exp(-xivec(jj)*time) - exp(muvec(ii)*time))/(xivec(jj) + muvec(ii)) ... 
                	     + (exp(xivec(jj)*time) - exp(muvec(ii)*time))/(xivec(jj) - muvec(ii))) ...
                   		/xivec(jj)/4;
                 
      		 Hprime(ii,jj,tt) = 0; 
      		 for kk =1:N
          		 T1 = - ((exp(-muvec(jj)*time) - exp(muvec(ii)*time))/(muvec(jj) + muvec(ii)) ...
                	        + (exp(xivec(kk)*time) - exp(muvec(ii)*time))/(xivec(kk) - muvec(ii))...
                 	          )/(muvec(jj)^2 - xivec(kk)^2)/2; 
              
           		 T2 = ((exp(xivec(kk)*time) - exp(muvec(ii)*time))/(xivec(kk) - muvec(ii)) ...
              		    + (exp(-xivec(kk)*time) - exp(muvec(ii)*time))/(xivec(kk) + muvec(ii)) ...
                               )/(muvec(jj) - xivec(kk))/4/xivec(kk); 
          
                 Hprime(ii,jj,tt) = Hprime(ii,jj,tt) + C_mat(ii,kk)*C_mat(jj,kk)*(T1 + T2); 
             end
	    end
	end
end

lam_mat = zeros(N,nt+1); 
sig_mat = zeros(N,nt+1); 
x_mat = zeros(N,nt+1); 
xp_mat = zeros(N,nt+1); 

lam_mat(:,1) = lam0; 
sig_mat(:,1) = sig0; 

for it =1:nt
    sig_mat(:,it+1) = diag(exp(-D2*tarray(it+1)))*sig0; 
    lam_mat(:,it+1) = diag(exp(-D1*tarray(it+1)))*lam0 + Tlam(:,:,it+1)'*sig0; 
    x_mat(:,it+1) = -G(:,:,it+1)*lam0 + H(:,:,it+1)*sig0; 
    xp_mat(:,it+1) = - Gprime(:,:,it+1)*lam0 + Hprime(:,:,it+1)*sig0; 
end


for it = 1:nt+1
    sig_comp_mat(:,it) = V2*sig_mat(:,it); 
    lam_comp_mat(:,it) = V1*lam_mat(:,it); 
    x_comp(:,it) = V1*x_mat(:,it); 
    xp_comp(:,it) = V2*xp_mat(:,it); 
end

costate_mat = [lam_comp_mat; sig_comp_mat]; 
state_mat = [x_comp; xp_comp]; 
% make plots to compare ---------------------------------------------------
figure(1)
subplot(1,2,1)
imagesc(tarray,1:2*N,round(state_mat - x',15)); colorbar; axis square
subplot(1,2,2)
imagesc(tarray,1:2*N,round(costate_mat - v(2*N+1:end,:),15)); colorbar; axis square

max(max(abs(costate_mat - v(2*N+1:end,:)))) 
max(max(abs(state_mat - x')))







