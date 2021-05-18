%  Main code to make plots on HM-nets 
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
kappa_arr = zeros(size(f_vec,2), size(g_vec,2)); 
L1_av_ctrl_arr = zeros(size(f_vec,2), size(g_vec,2)); 
L2_av_ctrl_arr = zeros(size(f_vec,2), size(g_vec,2)); 
L1_max_ctrl_arr = zeros(size(f_vec,2), size(g_vec,2)); 
L2_max_ctrl_arr = zeros(size(f_vec,2), size(g_vec,2)); 
kappa = 0; 

E1_sample_av = zeros(N,ntrial);     % array to collect E1 
E2_sample_av = zeros(N,ntrial);     % array to collect E1 
eig1_sample_av = zeros(N, ntrial);  % array to collect eigen values for layer 1
eig2_sample_av = zeros(N,ntrial);   % array to collect eigen values for layer 2

% degree heterogeneity and Cluster coeff, averaged over all trials
H_aver = zeros(size(f_vec,2), size(g_vec,2)); 
CC_aver = zeros(size(f_vec,2), size(g_vec,2)); 


for jfrac = 1:size(f_vec,2)
for jgam = 1:size(g_vec,2)
    kappa = 0;
for jtrial =1:ntrial
    [jtrial jfrac jgam]
    
    dirname = ['trial=',num2str(jtrial)];
    load(fullfile(dirname,['HMControl_ifrac=',num2str(jfrac),'_igam=',...
                            num2str(jgam),'itrial=',num2str(jtrial),'.mat']));
    
    kappa = kappa + dscale; 
    E1_sample_av(:,jtrial) = E1; 
    E2_sample_av(:,jtrial) = E2; 
    
    eig1_sample_av(:, jtrial) = xivec; 
    eig2_sample_av(:, jtrial) = muvec; 
      
    H_aver(jfrac, jgam) = H_aver(jfrac, jgam) + deg_het(jfrac,jgam); 
    CC_aver(jfrac, jgam) = CC_aver(jfrac, jgam) + clust_coeff(jfrac,jgam); 
    
end

kappa_arr(jfrac, jgam) = kappa/ntrial;

% all values of energy and eigenvalues averaged across all sample
% plot these quantities to get all values of 
aver_E1 = mean(E1_sample_av,2); 
aver_E2 = mean(E2_sample_av,2); 
aver_eig1 = mean(eig1_sample_av,2); 
aver_eig2 = mean(eig1_sample_av,2); 

% average and maximal control energies
L1_av_ctrl_arr(jfrac,jgam) = mean(aver_E1); 
L2_av_ctrl_arr(jfrac,jgam) = mean(aver_E2); 

L1_max_ctrl_arr(jfrac,jgam) = max(aver_E1); 
L2_max_ctrl_arr(jfrac,jgam) = max(aver_E2); 


end


end

H_aver = H_aver/ntrial; 
CC_aver = CC_aver/ntrial; 

figure(1)
imagesc(g_vec, f_vec, H_aver); 
axis xy; axis square
colormap('pink')
xlabel('\gamma');    ylabel('f')
set(gca,'fontsize',16)
title('Degree Heterogeneity', 'fontweight','normal')
colorbar

figure(2)
imagesc(g_vec, f_vec, CC_aver); 
axis xy
axis square
colormap('pink')
xlabel('\gamma');  ylabel('f')
set(gca,'fontsize',16)
colorbar
title('Average Clustering', 'fontweight','normal')

figure(3)
imagesc(g_vec, f_vec, L2_av_ctrl_arr); 
axis xy
axis square
colormap('pink')
xlabel('\gamma');  ylabel('f')
set(gca,'fontsize',16)
colorbar
title('Average control energy', 'fontweight','normal')

figure(4)
imagesc(g_vec, f_vec, L2_max_ctrl_arr); 
axis xy
axis square
colormap('pink')
xlabel('\gamma');   ylabel('f')
set(gca,'fontsize',16)
colorbar
title('Maximum control energy', 'fontweight','normal')










