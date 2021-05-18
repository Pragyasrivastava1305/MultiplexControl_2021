% author:       @pragyasr
% created in:    Summer, 2019
% modified in:   Spring, 2021

clear all;
close all;

% import colors
color_list; 
col_one = col_list(3,:);
fig_prop = 1;

% import python module
py.importlib.import_module('networkx');
N = 100;
T = Inf;
p = 0.1; 
directed = 0; 
ntrial = 40;

SymAdj = zeros(N,N);
int_vec = 2:2:N-1; 
real_vec = linspace(0.1,1,size(int_vec,2)); 

loop_vec = zeros(size(int_vec,2),4); 
loop_vec(:,1) = real_vec; 
loop_vec(:,2) = int_vec; 
loop_vec(:,3) = int_vec; 
loop_vec(:,4) = real_vec;
eig_array = zeros(N, ntrial); 

% to calculate distribution 
nbins  = 21; 

inet =2;

for inet =1:4
    
E_hist = zeros(nbins-1,1); 
xax = zeros(nbins-1,1);

for itrial = 1:ntrial
   
        if inet == 1
            G = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N)...
                                             ,0.25,pyargs('seed',int16(10*itrial),'directed',false)));
        elseif inet == 2
            G =net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(25),0.1)); 

        elseif inet == 3
            G = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(int16(N),int16(14)));

        elseif inet == 4
            G = net.helper.py_graph2adjmat(py.networkx.random_geometric_graph(int16(N),0.325)); 
        end
        G = G + G'; 
        
        % spectral properties
        [vv,ee] = eig(G); 
        
        eig_val = diag(ee)/ee(end,end);
        
        eig_array(:,itrial) = eig_val; 
        
        bin_edge = linspace(min(eig_val), max(eig_val),nbins); 
         
        h = histogram(eig_val,bin_edge,'Visible','off');
        hist_vec = h.Values; 
        hist_vec = hist_vec/sum(hist_vec); 
        
        E_hist = E_hist + hist_vec'; 
        
        med_edge = movmean(bin_edge,1); 
        
        xax = xax + med_edge(2:end)'; 
        
        % topological properties
        avg_clust(itrial , inet) = avgClusteringCoefficient(G); 
        deg_vec = sum(G); 
        
        deg_het(itrial , inet) = std(deg_vec); 
        chris_def(itrial, inet) = sum(sum(abs(deg_vec'...
                                            - deg_vec)))/(N*(N-1))/mean(deg_vec);
                                        
       first_moment( itrial, inet) = mean( eig_val ); 
       second_moment( itrial, inet) = var( eig_val ); 
       third_moment( itrial, inet) = skewness( eig_val ); 
       fourth_moment(itrial, inet) = kurtosis( eig_val );                                  
        
        
end
    
    E_hist = E_hist/ntrial; 
    
    % define bulk versus extremal values
    bulk_to_ex(inet) = sum(E_hist(11:20))/sum(E_hist(1:10)); 
    
    mean_Eval  = mean(eig_array');
    
    % spectrum
    mean_first_moment(inet) = mean( mean_Eval(:,inet) ); 
    mean_second_moment(inet) = var( mean_Eval(:,inet) ); 
    mean_third_moment(inet) = skewness( mean_Eval(:,inet) ); 
    mean_fourth_moment(inet) = kurtosis( mean_Eval(:,inet) ); 
    
    % topology
    mean_deg_het(inet) = mean(deg_het(:,inet)); 
    mean_chris_def(inet) = mean(chris_def(:,inet)); 
    mean_avg_clust(inet) = mean(avg_clust(:,inet)); 
    

    
    figure(inet)
    subplot(1,2,1)
    imagesc(G)
    set(gca,'fontsize',14);  colorbar
    xlabel('nodes'); ylabel('nodes')
    axis square
    
    
    subplot(1,2,2)
    bar(xax/ntrial , E_hist/ntrial , 'linewidth',2,'FaceAlpha',0.25);  hold all
    set(gca,'fontsize',14)
    ylim([ 0 1])
    title( ['extremal to bulk ratio = ', num2str(bulk_to_ex(inet))],'FontWeight','normal' )
   
    ylabel('P(\lambda)');  xlabel('\lambda')
    axis square; axis xy
    
                
    drawnow

%     saveas(gcf,['EigDist_inet=',num2str(inet),'.pdf'])
end


% fig for spectral and topological properties
if fig_prop ==1
    figure;
    
    subplot(3,2,1)
    violin(first_moment)
    title('mean eigenvalue')
    set(gca,'xticklabel',{'ER', 'WS', 'BA','RG'},'fontsize',10)
    
    subplot(3,2,2)
    violin(second_moment)
    title('spread of eigenvalue dist')
    set(gca,'xticklabel',{'ER', 'WS', 'BA','RG'},'fontsize',10)
    
    subplot(3,2,3)
    violin(third_moment) 
    title('skewness of eigenvalue dist')
    set(gca,'xticklabel',{'ER', 'WS', 'BA','RG'},'fontsize',10)
    
    
    subplot(3,2,4)
    violin(fourth_moment)
    title('kurtosis of eigenvalue dist')
    set(gca,'xticklabel',{'ER', 'WS', 'BA','RG'},'fontsize',10)
    
    subplot(3,2,5)
    violin(chris_def) 
    title('degree heterogeneity')
    set(gca,'xticklabel',{'ER', 'WS', 'BA','RG'},'fontsize',10)
    
    subplot(3,2,6)
    violin(avg_clust)
    title('average clustering')
    set(gca,'xticklabel',{'ER', 'WS', 'BA','RG'},'fontsize',10)
     
end














