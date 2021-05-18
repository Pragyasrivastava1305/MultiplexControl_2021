function [M] = init_duplex_identical_layers(inet1,rho1,N,opt_plot)
% function to create initial duplex
% INPUTS
% inet1: \in [1,4] to specify first layer topology
% rho1:  density of first layer
% N: size of the network (number of nodes)
% opt_plot: plot M if opt_plot=1. 
% OUTPUT: 
% M: duplex network

M  = zeros(2*N, 2*N); 
%  get generative paremeter for specified density
p1 = get_dens_param(rho1, inet1);  
% p2 = get_dens_param(rho2, inet2);  

% create input layer
if inet1 ==1
        G1 = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N),p1,pyargs('directed',false))); 
elseif inet1 ==2 
        G1 = net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(p1),0.1)); 
elseif inet1 ==3			 
    	G1 = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(int16(N),int16(p1)));
elseif inet1 == 4
    	G1 = net.helper.py_graph2adjmat(py.networkx.random_geometric_graph(int16(N),p1)); 
end


% make G1/G2 symmetric and traceless
G1 = (G1 + G1') - G1.*eye(N);
G2 = G1;
    
%create duplex
M(1:N,1:N) = G1;
M(N+1:2*N, 1+N:2*N) = G2;
M(N+1:2*N,1:N) = eye(N);

if opt_plot ==1
    imagesc(M)
    colorbar
end