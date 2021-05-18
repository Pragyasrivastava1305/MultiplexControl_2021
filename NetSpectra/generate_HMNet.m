function G = generate_HMNet(N, N_block, E, frac_within, gamma)
% Inputs: Total number of nodes N, number of nodes within each block
% N_block, total number of edges E, the fraction of edges frac_within to 
% add within blocks, and the scale-free exponent gamma.
%
% Ouput: NxN unweighted, undirected adjacency matrix G representing a
% hierarchically modular network.
%
% NOTE: We generate this network by combining the methods from
% 'generate_BlockNet_conn.m' and 'generate_StaticNet.m'. WE USE THIS
% FUNCTON TO GENERATE HM NETS FOR HUMAN INFO PROCESSING PAPER.

G = zeros(N);

alpha = 1/(gamma - 1);

num_blocks = floor(N/N_block);
Ns_block = [N_block*ones(1, num_blocks - 1), N_block + rem(N, N_block)];

% Connect the first node in block i to the first nodes in block i+1:
% num_blocks = N/N_block;
nodes_added = [1];

for i = 1:(num_blocks - 1)
    I = (i-1)*N_block + 1;
    J = i*N_block + 1;
    
    G(I,J) = 1;
    G(J,I) = 1;
    
    nodes_added = [nodes_added, J];
end

% Create random spanning tree starting at the first node in each block:
% E_tree = round(N_block*frac_within);
for i = 1:num_blocks
    E_tree = min([round(Ns_block(i)*frac_within), Ns_block(i) - 1]);
    
    for j = 2:(E_tree + 1)
        
        J = (i-1)*N_block + j;
        I = randi([(i-1)*N_block + 1, (i-1)*N_block + j - 1]);
        
        G(I,J) = 1;
        G(J,I) = 1;
        
        nodes_added = [nodes_added, J];
    end
end

% Add the remaining nodes in each block. Connect each to one of the nodes
% that's already been connected.
for i = 1:num_blocks
    N_block_temp = Ns_block(i);
    E_tree = round(N_block_temp*frac_within);
    
    for j = (E_tree+2):N_block_temp
        
        nodes_possible = setdiff(nodes_added, (i-1)*N_block + (1:N_block_temp));
        J = (i-1)*N_block + j;
        I = datasample(nodes_possible,1);
        
        G(I,J) = 1;
        G(J,I) = 1;
        
        nodes_added = [nodes_added, J];
    end
end

% Write down weights for each node:
W = (1:N).^(-alpha);
W_mat = W'*W;

% List the possible within- and between-community edges:
G_within = zeros(N);

for i = 1:num_blocks
    inds = (i-1)*N_block + (1:Ns_block(i));
    G_within(inds, inds) = 1;
end

edges_possible_within = find(triu(G_within - G > 0, 1));
edges_possible_between = find(triu(ones(N) - G_within - G > 0, 1));

% Select within- and between-community edges:

edge_weights_within = W_mat(edges_possible_within);
edge_weights_between = W_mat(edges_possible_between);

num_edges = E - N + 1;
% num_edges_within = round(frac_within*num_edges);
num_edges_within = round(frac_within*E - sum(sum(G_within.*G))/2);
if num_edges_within > 0
    edges_within = datasample(edges_possible_within, num_edges_within,...
        'Replace', false, 'Weights', edge_weights_within);
    
    G(edges_within) = 1;
end

num_edges_between = num_edges - num_edges_within;
if num_edges_between > 0
    edges_between = datasample(edges_possible_between, num_edges_between,...
        'Replace', false, 'Weights', edge_weights_between);
    
    G(edges_between) = 1;
end

G = double(G + G' > 0);

