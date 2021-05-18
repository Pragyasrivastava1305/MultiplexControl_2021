% This code matches solutions for state and costate variables.
% Contains comparisons at two levels: initial costate conditions (i)obtained
% numerically, and (ii) solved for numerically. 
% Variables are looped over rather than being in vectorised form. 
% pragyasr : Spring 2019

clear all;
close all;
% import python module
py.importlib.import_module('networkx');

% Dynamic parameters ------------------------------------------------------ 
% total time
T = 1; 
nt = 100; 
tarray =  linspace(0,T,nt+1); 
delt = T/nt ; 
Niter = 100; 

% create duplex by calling the following script. Ensure that the two 
% layers are not same
N = 50; 
CreateDuplex; 
muvec = D2;           xivec = D1; 

% ------------------- ALLOCATIONS/INITIALIZATION --------------------------
% initial condition for state variables
X1_init = zeros(N,1);    X2_init = zeros(N,1);     X0 = [X1_init; X2_init];

% input matrix
B = zeros(2*N,N);       B(1:N,1:N) = eye(N); 

% generate final state
iter =1; 
rng(iter);     X1_final = randn(N,1);
rng(2*iter);   X2_final = randn(N,1); 
XF = [X1_final; X2_final];              XF = XF/norm(XF);


% numerical solution to compare with -------------------------------------
[x,u,v,n_err] = min_eng_cont(Duplex, T, B, X0, XF, nt,0); 
costate_num = v(2*N+1:end,:); 
   
%% LEVEL 1 matching: extract initial condition from the numerical
% solution. Evolve using the analytical solution and match.

% Analytical solutions  -----------------------------------------------
L0_num = v(2*N+1:end,1);
lam0 = pinv(V1)*L0_num(1:N); 
sig0 = pinv(V2)*L0_num(N+1:end); 

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
% LEVEL 1 check complete [lam_comp, sig_comp] = v(2*N+1:end,:)


%% LEVEL 2 checks: extract initial condition of costate variable from
% numerical solution. Evolve state and costate variables using matrices of
% analytical solution

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

% calcualte control input and control energy 
ener_densL1 = zeros(nt+1,1); 
optim_u = -lam_comp_mat/2; 

for iter = 1:size(optim_u,2)
    ener_densL1(iter) = optim_u(:,iter)'*optim_u(:,iter); 
end

optim_EL1 = trapz(tarray,ener_densL1);

% LEVEL 2 checks complete ------------------------------------------------

%% LEVEL 3 check: Solve for initial value of costate variables and evolve 
% the initial conditions according to analytical solution

[state_mat2,costate_mat2] = AllModeSol_funV0(layer1,layer2,dscale,T,nt,XF); 
lam_comp_mat2 = costate_mat2(1:N,:); 

% calculate control input and control energy 
ener_densL2 = zeros(nt+1,1); 
optim_u = -lam_comp_mat2/2; 

for iter = 1:size(optim_u,2)
    ener_densL1(iter) = optim_u(:,iter)'*optim_u(:,iter); 
end

optim_EL2 = trapz(tarray,ener_densL1);


% make plots to compare ---------------------------------------------------
figure(2)
subplot(1,2,1)
imagesc(tarray,1:2*N,round(state_mat2 - x',15)); colorbar; axis square
subplot(1,2,2)
imagesc(tarray,1:2*N,round(costate_mat2 - v(2*N+1:end,:),15)); colorbar; axis square

max(max(abs(costate_mat2 - v(2*N+1:end,:)))) 
max(max(abs(state_mat2 - x')))


% LEVEL 2 gives better match ofcourse (10e-11) but Level 3 match is calculated
% completely analytically (10e-6 to 10e-9). 

