function [E1, E2, megaU_L1, megaU_L2, state_chk] = layer_control(Duplex, T, B, nt, V1, V2)
% Input
% Duplex:  Supra-adjacency matrix for duplex
% T: time horizon 
% B: Input matrix 
% nstep:  number of steps 
% V1, V2:  matrix with eigenvectors as columns for the two layers

% Output
% E1: optimal energies for layer 1 eigenvectors
% E2: same as E1, but for layer 2
% megaU_L1: a mega arrray that contains complete optimal input u 
%           megaU_L1: N x N*(nt+1) array
%           megaU_L1(i, :):  reshape(u_i, [N*(nt+1), 1])
%           here u_i is optimal input for i-th eigenvector, which has N x
%           (nt + 1) shap

% This code does not output the state variables but can be eaasily modified
% to do so. 

% Whether target state is being reached or not is saved in 'state_chk'
% vector. A non-zero entry of state_chk  implies that the system failed
% to reach the corresponding target eigenveector. This can happen when there are 
% numerical instabilities or when matrices are unstable (large positive
% eigenvalues). 

% megaU_L2: same as megaU_L1 but for layer 2
% state_chk: error between target eigenvector and the reached state    

tarray = linspace(0,T,nt+1);  

%  initial states
N = length(V1); 
X0 = zeros(2*N,1);
E1 = zeros(N,1); 
E2 = zeros(N,1);

megaU_L1 = zeros( N, N*(nt+1) ); 
megaU_L2 = zeros( N, N*(nt+1) );
state_chk = zeros(2*N, 1); 

for j_eig = 1:N
    
    % construct final state from the eigenvectors of layer 1
    XF1 = [V1(:,j_eig); zeros(N,1)];
    XF2 = [zeros(N,1); V2(:,j_eig)]; 
    
    % calculate state, costate and optimal input
    [x1,u1,v1,n_err1] = min_eng_cont(Duplex, T, B, X0, XF1, nt,0); 
    optim_u1 = u1(:,1:N)'; 

    [x2,u2,v2,n_err2] = min_eng_cont(Duplex, T, B, X0, XF2, nt,0); 
    optim_u2 = u2(:,1:N)'; 

    fin_E_chk1(j_eig) = round(max(abs(x1(end,:) - XF1')),10); 
    fin_E_chk2(j_eig) = round(max(abs(x2(end,:) - XF2')),10); 
            
    for iter = 1:nt+1
       ener_densL1(iter) = optim_u1(:,iter)'*optim_u1(:,iter); 
       ener_densL2(iter) = optim_u2(:,iter)'*optim_u2(:,iter); 
    end

    E1(j_eig) = trapz(tarray,ener_densL1);
    E2(j_eig) = trapz(tarray,ener_densL2);
            
    u1_vec = reshape(optim_u1,[N*(nt+1),1]); 
    u2_vec = reshape(optim_u2,[N*(nt+1),1]); 
    
    megaU_L1(j_eig,:) = u1_vec; 
    megaU_L2(j_eig,:) = u2_vec; 
            
end    % end of loop over eigenstates
    
 
state_chk(1:N) = fin_E_chk1; 
state_chk( N+1: 2*N) = fin_E_chk2; 
 

end   
    
    
    
    
    
    
    
    
    
    
    
    