%  Code to produce the average control energy vs kappa plots 
%  density is fixed, topology is fixed (the pair; WS-RG)

% author:       @pragyasr
% created in:    Summer, 2020

% final version of code which produces separate files for each trial
% Each trial result is saved in a different directory
% output files: 
% Topo_density_inet1='i'_inet2='j'_trial='itrial'.mat : saves mat file

% SupA_inet1='i'_inet2='j'_trial='itrial',_dens='idens'.csv : saves the network 

% optimU_inet1='i'_inet2='j'_trial='itrial',_dens='idens'.csv : saves optim_u for layer 1


clear all;
% close all;

N = 100;
T = 5;
nt =50;
tarray = linspace(0,T,nt+1); 
p = 0.1; 
directed = 0; 
ntrial = 5; 
 
% assign control input matrix
B = zeros(2*N,N); 
B(1:N,1:N) = eye(N); 
kappa_Vec = 1:0.1:10; 

% vector to loop over to increase connectivity/density of second layer for WS/BA
int_vec = [2:10:N-1 N-1];  

% vector to loop over to increase connectivity/density of second layer ER/RG
real_vec = linspace(0.1,1,size(int_vec,2));  % parameter to create ER and RG nets

E1_array = zeros(N, size(kappa_Vec,2)); 
E2_array = zeros(N, size(kappa_Vec,2)); 
optimU_L1 = zeros(N, (nt+1)*N); 
optimU_L2 = zeros(N, (nt+1)*N); 

rho0 = 0.25;
inet1 = 2;
inet2 = 4;
idens = 6; 


for inet1 =1:4
for inet2=1:4
        
for itrial = 1:ntrial
    
for ik =1:size(kappa_Vec,2)

    kappa = kappa_Vec(ik);
    
    %for inet1 = 1:4    % first layer loop 
  
    p = get_dens_param(rho0, inet1); 

    % create first layer
    if inet1 ==1
            G1 = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N),p,pyargs('directed',false))); 
    elseif inet1 ==2 
            G1 = net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(p),0.1)); 
    elseif inet1 ==3			 
            G1 = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(int16(N),int16(p)));
    elseif inet1 == 4
            G1 = net.helper.py_graph2adjmat(py.networkx.random_geometric_graph(int16(N),0.325)); 
    end
    
    % make G1 symmetric and traceless
    G1 = (G1 + G1') - G1.*eye(N);
    
    % calculate the density of network
    rho1 = sum(sum(G1))/N/(N-1); 
    
    %create duplex
    Duplex0(1:N,1:N) = G1;
    Duplex0(N+1:2*N,1:N) = kappa*eye(N);
   

    %for inet2 = 1:4   % second layer loop 
    
    %mkdir(['trial=',num2str(itrial),'/inet1=',num2str(inet1),'_inet2=',num2str(inet2)])
    %dirname = ['trial=',num2str(itrial),'/inet1=',num2str(inet1),'_inet2=',num2str(inet2)]; 
    
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
    
    % initiate loop over the density of second layer network
    %for idens = 1:size(loop_vec,2)
   	[ itrial, kappa ]
	%  if idens ==2;
	%      break 
	%  end 
   
 	% create second layer with increasing density
    if inet2 ==1
        G2 = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N),loop_vec(idens),pyargs('directed',false))); 
    elseif inet2 ==2 
        G2 = net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(loop_vec(idens)),0.1)); 
    elseif inet2 ==3 
        G2 = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(int16(N),int16(loop_vec(idens))));
    elseif inet2 ==4
        G2 = net.helper.py_graph2adjmat(py.networkx.random_geometric_graph(int16(N),loop_vec(idens))); 
    end
      
       
     	G2 = (G2 + G2') -  G2.*eye(N);  
      	Duplex0(N+1:2*N, N+1:2*N) = G2;  
       
		% calculate density of second layer
      	rho2 = sum(sum(Duplex0(N+1:2*N,1+N:2*N)))/N/(N-1); 
   
  		% stabilize the network if required
        max_eig_G1 =  eigs(G1,1);
        
        % non-dimensionalize time by dividing by the largest eigenvalye of G1
        Duplex = Duplex0/max_eig_G1; 
        
       
		% extract layer networks
      	layer1 = Duplex(1:N,1:N); 
       	layer2 = Duplex(N+1:2*N, N+1:2*N); 
      	off_mat = Duplex(N+1:end,1:N); 

        % write duplex matrix to analyse later 
        %writematrix(Duplex, fullfile(dirname,['SupA_inet1=',num2str(inet1),'_inet2=',...
        %                                 num2str(inet2),'_trial=',num2str(itrial),'_dens=',num2str(idens),...
        %                                                                        '.csv']),'delimiter','tab');  
            
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
              
        % calculate the state, costate and energy corresponding to each eigendirection
        X0 = zeros(2*N,1); 
        
        for j_eig = 1:N
%           j_eig
            % specify final states along eigen-directions
            XF1 = [V1(:,j_eig); zeros(N,1)]; 
            XF2 = [zeros(N,1); V2(:,j_eig)]; 
           
            [x1,u1,v1,n_err1] = min_eng_cont(Duplex, T, B, X0, XF1, nt,0); 
            optim_u1 = u1(:,1:N)'; 

            [x2,u2,v2,n_err2] = min_eng_cont(Duplex, T, B, X0, XF2, nt,0); 
            optim_u2 = u2(:,1:N)'; 

            fin_E_chk1(j_eig,ik) = round(max(abs(x1(end,:) - XF1')),10); 
            fin_E_chk2(j_eig,ik) = round(max(abs(x2(end,:) - XF2')),10); 
            
            
            
            for iter = 1:nt+1
                ener_densL1(iter) = optim_u1(:,iter)'*optim_u1(:,iter); 
                ener_densL2(iter) = optim_u2(:,iter)'*optim_u2(:,iter); 
            end

            E1(j_eig) = trapz(tarray,ener_densL1);
            E2(j_eig) = trapz(tarray,ener_densL2);
            
%           % write control inputs to analyse later 
            u1_vec = reshape(optim_u1,[N*(nt+1),1]); 
            u2_vec = reshape(optim_u2,[N*(nt+1),1]); 
            
            optimU_L1(j_eig,:) = u1_vec; 
            optimU_L2(j_eig,:) = u2_vec; 
            
           
        end    % end of loop over eigenstates
        
        
        
%         if mod(itrial,5) ==0
%         writematrix(optimU_L1, fullfile(dirname,['optimU_L1_inet1=',num2str(inet1),...
%                                     '_inet2=',num2str(inet2),'_trial=',num2str(itrial),...
%                                     '_dens=',num2str(idens),'.csv']),'delimiter','tab') ;
%                                 
%         writematrix(optimU_L2, fullfile(dirname,['optimU_L2_inet1=',num2str(inet1),...
%                                     '_inet2=',num2str(inet2),'_trial=',num2str(itrial),...
%                                     '_dens=',num2str(idens),'.csv']),'delimiter','tab') ;
%         end
%         
        E1_array(:,ik) = E1;
        E2_array(:,ik) = E2;
        
    
 %end   % end of loop over density
 

%end   % end of inet2 loop
 
%end   % end of inet1 loop
 
%end   % end of trial loop

end
%save(fullfile(dirname,['topodensity_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'itrial=',num2str(itrial),'.mat']))
save(['topodensity_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'itrial=',num2str(itrial),'.mat'])

end 
end
end