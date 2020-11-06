% Code to create duplex in different topological classes and densities
% author:       @pragyasr
% created in:    Summer, 2020

clear all;
%close all;
% import python module
%py.importlib.import_module('networkx');

N = 100;
p = 0.1; 
directed = 0; 
norm = 0; 
nondim = 1;

% assign control input matrix
B = zeros(2*N,2*N); 
B(1:N,1:N) = eye(N); 

% vector to loop over to increase connectivity/density of second layer for WS/BA
int_vec = 2:4:N-1;  

% vector to loop over to increase connectivity/density of second layer ER/RG
real_vec = linspace(0.1,1,size(int_vec,2));  % parameter to create ER and RG nets

inet1 = 1; 

inet2 = 3; 

% for inet1 =1:4    % first layer loop
% for inet2 =1:4    % second layer loop 

    % specify the class of networks in first layer and second layer
    % 1 == ER, 2 == WS, 3 == BA, 4 == RG
  
    k2 = floor(N/4); % parameter to set density of WS in first layer
    k4 = floor(N/8); % parameter to set density of BA in first layer

    % create first layer
    if inet1 ==1
        G1 = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N),0.25,pyargs('directed',false))); 
    elseif inet1 ==2 
        G1 = net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(k2),0.1)); 
    elseif inet1 ==3			 
    	G1 = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(int16(N),int16(k4)));
    elseif inet1 == 4
    	G1 = net.helper.py_graph2adjmat(py.networkx.random_geometric_graph(int16(N),0.25)); 
    end

    % set loop_vec to a vetor which is used to create second layer of varying density
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
    % create duplex
    Duplex0(1:N,1:N) = G1;
    Duplex0(N+1:2*N,1:N) = eye(N);
    idens = 19; 
	
    % initiate loop to change density of second layer network
  
   
    	sd = randperm(100000,1);  % random seed
    
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
       
    rho2 = sum(sum(Duplex0(N+1:2*N,1+N:2*N)))/N/(N-1); 
        
	Duplex = Duplex0; 
	% stabilize if norm =1
    if norm ==1 
     	max_eig = max(real(eig(Duplex0))); 
  		Duplex = Duplex0/(max_eig +1) -1 ; 
    end

    % non-dimensionalize time by dividing by the largest eigenvalye of G1 
	% if nondim = 1
    if nondim ==1 
		Duplex = Duplex0/eigs(G1,1); 
    end  

%end
%end




