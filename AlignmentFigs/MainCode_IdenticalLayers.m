% code to calculate optimal energies for duplexes with identical layers
clear all; close all

opt_mov =0; 
opt_plt =0; 
opt_norm =0; 
py.importlib.import_module('networkx');

N = 100;     % network size
m = 1;       % number of modes constituting the target state
T = 5;       % time horizon 
nt = 50;     % number of time points
tarray =  linspace(0,T,nt+1);    % time array
Ntrial = 100;                    % number of trials for each duplex
B = zeros(2*N, N);               % input matrix
B(1:N,1:N) =  eye(N); 

m2 = zeros(m,4); 

if opt_mov ==1 
    mov = VideoWriter('Test_movie.avi');  % initiate movie 
end

% create original duplex by calling the following script
rho1 = 0.25; 
rho2 = rho1;
% col = [184, 115, 51]/255;
color = [0.8, 0.8, 0.8]; 

inet1 = 1; 
for inet1 =1:4
 
     if inet1 ==1
%       col_sch = col1;
        col_one = [95 158 160]/255;
    elseif inet1 ==2
%       col_sch = col2;
        col_one = [186,181,147]/255;
    elseif inet1 == 3
%       col_sch = col3;
        col_one = [167,196,139]/255;
    else
%       col_sch = col4;
        col_one = [231,130,162]/255; 
    end
    
    
    Duplex0 = init_duplex_identical_layers(inet1,rho1,N, 0); 

    scale_fac = eigs(Duplex0(1:N, 1:N),1 );   % scale by largest eigenvalue of first layer

    Duplex = Duplex0/scale_fac; 

    L1 = Duplex(1:N, 1:N); 
    % set the two layers to be the same
    L2 = Duplex(1+N:2*N, 1+N:2*N); 

    % eigenvalues and eigenmodes
    [V1,D1] = eig(L1);      
    [D1,I1] = sort(diag(D1));  V1 = V1(:,I1); 
    P = V1; 
    Q = P;
    xi = D1; 
    mu = xi; 

    % Set initial conditions
    X1_init = zeros(N,1);    
    X2_init = zeros(N,1);
    X0 = [X1_init; X2_init];

    omega_21 = zeros(nt+1,N); 
    E_21 = zeros(nt+1,N); 
 
    % pick the eigenvectors that contribute to final state @random
    m2(:,inet1) = sort(randi(N,[1,m]));
    
    vec_P = zeros(N,1);     
    vec_Q = zeros(N,1);
    rvec = 0.5*rand(m,1) + 0.5;

    % assign the eigenvectors of second layer as the target-modes
    for col =1:N
        col
        vec_Q = Q(:,col); 
        vec_Q = vec_Q/norm(vec_Q); % normalize the final state
        XF = [vec_P; vec_Q]; 
        XF = XF/norm(XF);

        [state,u,v,n_err] = min_eng_cont(Duplex, T, B, X0, XF, nt,0); 

        % calculate omega_21
        for iter = 1:nt+1
            omega_21(iter,:) = u(iter,:)*P; 
            
            % element-wise square of omega to calculate energy in channels
            E_21(iter,:) = omega_21(iter,:).^2;
            
            % energy at each time in original and projected spaces (they should be same)
            dens_proj_space(iter) = omega_21(iter,:)*omega_21(iter,:)'; 
            ener_dens(iter) = u(iter,:)*u(iter,:)';  
        end
    
        % these two quantities should be same: it checks out
        E_tot(col) = trapz(tarray,ener_dens); 
        E_chk(col) = trapz(tarray,dens_proj_space); 
        
       
       
        % identify the channel that carries the most energy and the
        % fraction of energy. For this sum over the time dimension
        [E21_max(col), ind_max(col)] = max(sum(E_21)); 
        E_frac(col) = E21_max(col)/sum(sum(E_21));
        
        
 end
   
subplot(1,3,1)
s1 = scatter(mu,E_tot);  hold on
s1.MarkerEdgeColor = 'k'; 
s1.MarkerFaceColor = col_one;
s1.SizeData = 60;
set(gca,'fontsize',14)
box on
axis square

subplot(1,3,2)
s2=scatter(1:N, ind_max); hold on 
s2.MarkerEdgeColor = 'k'; 
s2.MarkerFaceColor = col_one;
set(gca,'fontsize',14)
s2.SizeData = 60; 

box on
axis square

subplot(1,3,3)
s3 = scatter(1:N, E_frac); hold on 
s3.MarkerEdgeColor = 'k'; 
s3.MarkerFaceColor = col_one;
set(gca,'fontsize',14)
s3.SizeData = 60; 

box on
axis square

drawnow

end













