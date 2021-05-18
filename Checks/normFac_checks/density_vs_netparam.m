% author:       @pragyasr
% created in:    Summer, 2019
clear all;
%close all;
% import python module
py.importlib.import_module('networkx');
N = 100;
T = Inf;
p = 0.1; 
directed = 0; 
ntrial = 20;

ER1 = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N),0.1,pyargs('directed',false))); 
ER1 = 0.5*(ER1 + ER1'); 


SymAdj = zeros(N,N);
int_vec = 2:2:N-1; 
real_vec = linspace(0.1,1,size(int_vec,2)); 

loop_vec = zeros(size(int_vec,2),4); 
loop_vec(:,1) = real_vec; 
loop_vec(:,2) = int_vec; 
loop_vec(:,3) = int_vec; 
loop_vec(:,4) = real_vec;

for inet =1:4

%     if inet ==1
%         net_opt = 'ER';
%         loop_vec = real_vec; 
%     elseif inet ==2 
%         net_opt = 'WS';
%         loop_vec = int_vec; 
%     elseif inet ==3
%         net_opt = 'BA';
%         loop_vec = int_vec; 
%     else 
%         net_opt = 'RG';
%         loop_vec = real_vec; 
%     end
        
for idens = 1:size(loop_vec,1)
    idens
    for itrial = 1:ntrial
   
        if inet == 1
            G = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N)...
                                             ,loop_vec(idens,inet),pyargs('seed',int16(idens),'directed',false)));
        elseif inet == 2
            G =net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(loop_vec(idens,inet)),0.1)); 

        elseif inet == 3
            G = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(int16(N),int16(loop_vec(idens,inet))));

        elseif inet == 4
            G = net.helper.py_graph2adjmat(py.networkx.random_geometric_graph(int16(N),loop_vec(idens,inet))); 
        end
    
    
        % Symmetrize
        SymAdj = G;     
        SymAdj = (SymAdj + SymAdj') - SymAdj.*eye(N); 

        %imagesc(SymAdj); drawnow
        % get degree of each node and the average degree of network
%         deg_vec = sum(SymAdj);
%         av_deg(idens,inet) = mean(deg_vec); 
%         max_deg(idens,inet) = max(deg_vec);
%         std_deg(idens,inet) = std(deg_vec); 

        % network density 
        rho_trial(itrial) =  sum(sum(SymAdj))/N/(N-1); 
        
    end
    
       rho(idens,inet) = mean(rho_trial); 
end

figure(inet)
plot(loop_vec(:,inet), rho(:,inet),'o-.','LineWidth',2); hold on
drawnow

end