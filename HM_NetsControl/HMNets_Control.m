%  Main code to calculate control properties of HM-nets
%  For each duplex network, 'ntrial' number of trials are generated, 
%  and Mat files are saved for each. Finally, supra-adjacency matrices do
%  not need to be saved separately, as they get saved in mat files. 

clear all;     close all;

% import colors
color_list; 
col_one = col_list(3,:);

% Set network parameters
N = 100;           % number of nodes
N_block = 10;      % number of N-block
E = 500;           % E is number of edges   
B = zeros(2*N,N); 
B(1:N,1:N) = eye(N); 

f_vec= linspace(0.1, 0.9, 5); 
g_vec = linspace(2.1,3, 5);

% simulation parameters
ntrial = 20; 
T = 5; 
rho0 = 0.25; 
nt =50; 
tarray = linspace(0,T,nt+1); 

% initilization of arrays 
eig_array = zeros(N,ntrial); 

for itrial =2:ntrial
      
    mkdir(['trial=',num2str(itrial)]);
    dirname = ['trial=',num2str(itrial)];
    
    p = get_dens_param(rho0, 1); % construct first layer frrom ER  
    
    G1 = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N),p,pyargs('directed',false)));
    G1 = (G1 + G1') - G1.*eye(N);
    
    % calculate the density of network
    rho1 = sum(sum(G1))/N/(N-1); 
    
    %create duplex
    Duplex0(1:N,1:N) = G1;
    Duplex0(N+1:2*N,1:N) = eye(N);
   
    for ifrac = 1:size(f_vec,2)
    for igam = 1:size(g_vec,2)
    
        [itrial ifrac igam]
    
        G2 = generate_HMNet(N, N_block, E, f_vec(ifrac), g_vec(igam)); 
        Duplex0(N+1:2*N, N+1: 2*N) = G2; 
    
        % topological properties 
        deg_vec = sum(G2);
        rho_vec(ifrac, igam) = sum(sum(G2))/N/(N-1); 
        clust_coeff(ifrac, igam) = avgClusteringCoefficient(G2); 
        deg_het(ifrac, igam) = sum(sum(abs(deg_vec'...
                        - deg_vec)))/(N*(N-1))/mean(deg_vec);
    
    
        max_eig_G1 =  eigs(G1,1);
        
        % non-dimensionalize time by dividing by the largest eigenvalye of G1
        Duplex = Duplex0/max_eig_G1; 
   
        % extract layers of scaled network
        layer1 = Duplex(1:N,1:N); 
        layer2 = Duplex(N+1:2*N, N+1:2*N); 
        off_mat = Duplex(N+1:end,1:N); 

        % extract scaling of off-diagonal connections  
        off_diag = diag(off_mat); 
        dscale = off_diag(1);  
  
        %------------------- eigen-decomposition of each layer --------------------
        [V1,D1] = eig(layer1);      
        [V2,D2] = eig(layer2); 
        [D1,I1] = sort(diag(D1));  V1 = V1(:,I1); 
        [D2,I2] = sort(diag(D2));  V2 = V2(:,I2); 

        xivec = D1; 
        muvec = D2; 
              
        % calculate the state, costate and energy corresponding to each
        % eigendirection: X0  is specified in function itself.
     
        [E1, E2, megaU_L1, megaU_L2, state_chk] = layer_control(Duplex, T, B, nt, V1, V2); 
     
        save(fullfile(dirname,['HMControl_ifrac=',num2str(ifrac),'_igam=',num2str(igam),'itrial=',num2str(itrial),'.mat'])) 
    
    end
%     figure; 
%     colormap('bone')
%     imagesc(Duplex); 
%     drawnow
 
    
    end
    
    
 end







