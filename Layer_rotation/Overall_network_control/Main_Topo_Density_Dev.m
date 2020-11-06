% Code to calculate control energies/trajectories for each layer of the
% duplex network, for varying densities and topologise. 
% author:       @pragyasr
% created in:    Summer, 2020

clear all;
% close all;

N = 100;
T = 5;
nt =50;
tarray = linspace(0,T,nt+1); 
p = 0.1; 
directed = 0; 
ntrial = 1; 
 
% assign control input matrix
B = zeros(2*N,2*N); 
B(1:N,1:N) = eye(N); 

% vector to loop over to increase connectivity/density of second layer for WS/BA
int_vec = 2:16:N-1;  

% vector to loop over to increase connectivity/density of second layer ER/RG
real_vec = linspace(0.1,1,size(int_vec,2));  % parameter to create ER and RG nets

xi_array = zeros(N*ntrial, size(int_vec,2)); 
mu_array = zeros(N*ntrial, size(int_vec,2)); 
E1_array = zeros(N*ntrial, size(int_vec,2)); 
E2_array = zeros(N*ntrial, size(int_vec,2)); 

clen = 9; 
blue = [0,0,1];
green = [0, 1, 0];
colors_p = [linspace(blue(1),blue(1),clen)', linspace(blue(2),green(2),clen)', linspace(blue(3),green(3),clen)'];

% dirname = 'WS_eigen_analysis'; 

dirname = 'WS_randomization'; 
rho0 = 0.25;
inet1 =2;

% for inet1 = 1:4    % first layer loop

p_rewire_vec = linspace(0.2,1,5); 

for j_rew = 1:size(p_rewire_vec,2)
    
    p_rew = p_rewire_vec(j_rew); 
    
for inet2 = 1:4    % second layer loop 
 

    % specify the class of networks in first layer and second layer
    % 1 == ER, 2 == WS, 3 == BA, 4 == RG
  
%     k2 = floor(N/4); % parameter to set the density of WS in first layer
%     k4 = floor(N/8); % parameter to set the density of BA in first layer

    p = get_dens_param(rho0, inet1); 

    % create first layer
    if inet1 ==1
        G1 = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N),p,pyargs('directed',false))); 
    elseif inet1 ==2 
        G1 = net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(p),p_rew)); 
    elseif inet1 ==3			 
    	G1 = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(int16(N),int16(p)));
    elseif inet1 == 4
    	G1 = net.helper.py_graph2adjmat(py.networkx.random_geometric_graph(int16(N),p)); 
    end

    % select loop_vec to create second layer of varying density
    if inet2 ==1
    		loop_vec = real_vec; 
    elseif inet2 ==2
   		 loop_vec = int_vec; 
    elseif inet2 ==3
   		 loop_vec = int_vec; 
    elseif inet2 ==4
    		loop_vec = real_vec; 
    end

    % make G1 symmetric and traceless
    G1 = (G1 + G1') - G1.*eye(N);
    
    % calculate the density of network
    rho1 = sum(sum(G1))/N/(N-1); 
    
    %create duplex
    Duplex0(1:N,1:N) = G1;
    Duplex0(N+1:2*N,1:N) = eye(N);
    %idens = 4; 
    
    % initiate loop over the density of second layer network
    for idens = 1:size(loop_vec,2)
   	[inet1 inet2 idens] 
	%  if idens ==2;
	%      break 
	%  end 
   
    
	% generate several iteration of second layer
	for itrial = 1:ntrial
        if itrial ==2
            break
        end
        
    	sd = randperm(500000,1);  % random seed
    
   		% create second layer with increasing density
       	if inet2 ==1
           		G2 = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N),loop_vec(idens),pyargs('directed',false))); 
        elseif inet2 ==2 
                G2 = net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(loop_vec(idens)),0.1,int16(sd))); 
      	elseif inet2 ==3 
           		 G2 = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(int16(N),int16(loop_vec(idens)),int16(sd)));
       	elseif inet2 ==4
         		 G2 = net.helper.py_graph2adjmat(py.networkx.random_geometric_graph(int16(N),loop_vec(idens))); 
        end
      
       
     	G2 = (G2 + G2') -  G2.*eye(N);  
      	Duplex0(N+1:2*N, N+1:2*N) = G2;  
       
		% calculate density of second layer
      	rho2(idens,itrial) = sum(sum(Duplex0(N+1:2*N,1+N:2*N)))/N/(N-1); 
   
  		% calculate spectral distance between two layers
      	dspec0(idens,itrial) = spec_dist(G1, G2, 1);
       
  		% stabilize the network if required
        
     	% max_eig = max(real(eig(Duplex0))); 
  		% Duplex = Duplex0/(max_eig +1) -1 ; 
        max_eig_G1(idens) =  eigs(G1,1);
        
        % non-dimensionalize time by dividing by the largest eigenvalye of G1
        Duplex = Duplex0/max_eig_G1(idens); 
          
 		% Choose this option for no scaling and no stabilization
      	% Duplex = Duplex0;  
        
		% extract layer networks
      	layer1 = Duplex(1:N,1:N); 
       	layer2 = Duplex(N+1:2*N, N+1:2*N); 
      	off_mat = Duplex(N+1:end,1:N); 

        % write duplex matrix to analyse later 
        writematrix(Duplex, fullfile(dirname,['Duplex_adj_inet1=',num2str(inet1),...
                                '_inet2=',num2str(inet2),'_idens=',num2str(idens),...
                                  '_itrial=',num2str(itrial), '.csv']),'delimiter','tab') ; 
        
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
             
        % arrays to save xi and mu values for all the trials and all
        % densities. The length of vector and ntrials are reshaped into
        % one-dim. To retrieve, reshape xi_array(:,idens) =
        % xi_array(N,ntrial,idens)
         
        xi_array((itrial-1)*N+1: itrial*N, idens) = xivec;
        mu_array( (itrial-1)*N+1: itrial*N, idens) = muvec;
                    
         
        % calculate the state, costate and energy corresponding to each eigendirection
        X0 = zeros(2*N,1); 
        
        for j_eig = 1:N
            
            %j_eig
            % specify final states along eigen-directions
            XF1 = [V1(:,j_eig); zeros(N,1)]; 
            XF2 = [zeros(N,1); V2(:,j_eig)]; 
           
            [x1,u1,v1,n_err1] = min_eng_cont(Duplex, T, B, X0, XF1, nt,0); 
            optim_u1 = u1(:,1:N)'; 

            [x2,u2,v2,n_err2] = min_eng_cont(Duplex, T, B, X0, XF2, nt,0); 
            optim_u2 = u2(:,1:N)'; 

            for iter = 1:nt+1
                ener_densL1(iter) = optim_u1(:,iter)'*optim_u1(:,iter); 
                ener_densL2(iter) = optim_u2(:,iter)'*optim_u2(:,iter); 
            end

            E1(j_eig,itrial) = trapz(tarray,ener_densL1);
            E2(j_eig,itrial) = trapz(tarray,ener_densL2);
            
            % write control inputs to analyse later 
            if mod(j_eig,10)==0
                new_mat = [u1(:,1:N), u2(:,1:N)]; 
                
                writematrix(new_mat, fullfile(dirname,['Duplex_adj_inet1=',num2str(inet1),...
                                '_inet2=',num2str(inet2),'_idens=',num2str(idens),...
                                  '_itrial=',num2str(itrial),'_jeig=',num2str(j_eig), '.csv']),'delimiter','tab') ; 
            end

        end    % end of loop over eigenstates
        
        E1_array((itrial-1)*N+1: itrial*N, idens) = E1(:,itrial);
        E2_array((itrial-1)*N+1: itrial*N, idens) = E2(:,itrial);
        
 end        % end of loop over trials of 2nd layer
    
    
    
end         % end of loop over density
 
 save(fullfile(dirname,['EigenEnergy_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'rew=',num2str(p_rew),'_original.mat']))

      
end   % end of inet2 loop
end   % end of inet1 loop

  
 
 
 
 
 
 
