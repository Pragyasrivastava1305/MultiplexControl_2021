function [state,costate] = AllModeSol_funV0(A1,A2,dscale,T,nt,XF) 
% This function calculates the state and costate trajectories for a duplex
% with linear dynamics in continuous time. The function solves for state
% and costate trajectories in a representation in which each layer of the
% duplex is diagonal. Then the solutions are transformed back to the
% original basis.  Origin is assumed to be initial condition. A1 and A2 are
% assumed to be extracted from an already stabilized Duplex.
% @pragya: Fall 2019

% INPUT : -----------------------------------------------------------------
% A1	: Adjacency matrix of layer 1
% A2	: Adjacency matric of layer 2
% dscale: strength of interlayer connections
% T 	: time horizon
% nt	: number of intervals in [0,T] interval

% OUTPUT : ----------------------------------------------------------------
% state  : 2N-by-1 dimensional state trajectories
% costate: 2N-by-1 dimensional costate trajectories
% -------------------------------------------------------------------------
% A1 = Sym_WS; 
% A2 = Sym_BA; 

tarray = linspace(0,T,nt+1); 
N = length(A1); 

% Get eigenvalues and alignment matrix ------------------------------------
[V1,D1] = eig(A1); 
[V2,D2] = eig(A2); 

% sort in decreasing order of magnitude (since all eigenvalues at 
% this state will be negative, the most negative eigenvalues will
% be the first)
[D1,I1] = sort(diag(D1)); 
[D2,I2] = sort(diag(D2));
V1 = V1(:,I1); 
V2 = V2(:,I2); 
xivec = D1; 
muvec = D2;
C_mat = zeros(N,N);

for i =1:N
    for j=1:N 
      C_mat(i,j) = V2(:,i)'*V1(:,j); 
    end
end

% scale C_mat down by apprpriate factor (obtained from the scaling factor by which 
% off-diagonal matrix elements are scaled down due to normalization and stabilization of
% original matrix)
C_mat = dscale*C_mat; 


% Define matrices for evolution -------------------------------------------
% Allocate matrices
Tlam = zeros(N,N,nt+1);
G = zeros(N,N,nt+1); 
H = zeros(N,N,nt+1); 
Gprime = zeros(N,N,nt+1); 
Hprime = zeros(N,N,nt+1); 

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

% -------------- obtain initial conditions of costate variables  ----------
X1F = pinv(V1)*XF(1:N); 
X2F = pinv(V2)*XF(N+1:end); 

% get all the matrices at terminal time
GT = G(:,:,end); 
HT = H(:,:,end); 
GprimeT = Gprime(:,:,end); 
HprimeT = Hprime(:,:,end); 

% determine initial conditions of costate variables
mat_int1 = HprimeT*pinv(HT); 

den_lam = -GprimeT + mat_int1*GT;

num_lam =  (X2F - mat_int1*X1F);

lam0 = pinv(den_lam)*num_lam; 

sig0 = pinv(HT)*( X1F + GT*lam0 ); 
size(sig0);
 
% Solve for state and costate variables -----------------------------------
% allocate intermediate state and costate variables 
state_int = zeros(2*N,nt+1); 
costate_int = zeros(2*N,nt+1); 
costate_int(1:N,1) = lam0; 
costate_int(N+1:end,1) = sig0; 

for it =1:nt
    costate_int(N+1:end, it+1) =  diag(exp(-D2*tarray(it+1)))*sig0; 
    
    costate_int(1:N,it+1) = diag(exp(-D1*tarray(it+1)))*lam0 + Tlam(:,:,it+1)'*sig0; 
    
    state_int(1:N,it+1) = - G(:,:,it+1)*lam0 + H(:,:,it+1)*sig0; 
    
    state_int(N+1:end,it+1) =  - Gprime(:,:,it+1)*lam0 + Hprime(:,:,it+1)*sig0;
   
end


% transform back to original basis
state = zeros(2*N,nt+1); 
costate = zeros(2*N,nt+1); 
for it = 1:nt+1
    state(1:N,it) = V1*state_int(1:N,it); 
    state(1+N:2*N,it) = V2*state_int(1+N:2*N,it); 
    costate(1:N,it) = V1*costate_int(1:N,it); 
    costate(N+1:2*N,it) = V2*costate_int(N+1:2*N,it); 
end

end


















