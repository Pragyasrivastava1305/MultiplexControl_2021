clear all
close all;
ndens =25; 
num_trial =40;
ysc = 2;  % if ysc =1, linear scale, if ysc =2, log scale.
N = 100;
ntop =4;
% code to define colors
def_colors; 

% plot options
Adj_plot =0;
Eig_plot =0;
fig_save =0;

fig_dir = 'Figures/Eigval_checks'; 

inet1 = 1; 
inet2 = 2;
G_old = zeros(2*N,2*N); 


for inet1 =1:ntop
 
for inet2  = 1:ntop
    
    rho2_array = zeros(ndens,num_trial); 
    EigL1_dens_arr = zeros(N,ndens); 
    EigL2_dens_arr = zeros(N,ndens); 
    
    if Eig_plot ==1
        
    end
    for jdens =1:25 
        
    EigL1_trial_arr = zeros(N, num_trial);
    EigL2_trial_arr = zeros(N,num_trial);
   
    for jtrial =1:40
        
        dirname = ['trial=',num2str(jtrial),'/inet1=',num2str(inet1),'_inet2=',num2str(inet2)]; 

        % load mat file
        load(fullfile(dirname,['topodensity_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'itrial=',num2str(jtrial),'.mat']))

        rho2_array(:,jtrial) =  rho2; 

    %   for idens = 1: %size(loop_vec,2)
        SupraA = load(fullfile(dirname,['SupA_inet1=',num2str(inet1),'_inet2=',...
                                                 num2str(inet2),'_trial=',num2str(jtrial),'_dens=',num2str(jdens),...
                                                                                  '.csv']));  

        % plot the difference in supraadjacency matrix between trials to check that 
        % each trial has a different random network
        if Adj_plot==1
                subplot(8,5,jtrial)
                imagesc(SupraA - G_old); drawnow
                G_old = SupraA;
        end

        if Eig_plot ==1
          % extract layers
          layer1 = SupraA(1:N,1:N); 
          layer2 = SupraA(N+1:2*N, N+1:2*N); 

          [VL1,DL1] = eig(layer1); 
          [VL2,DL2] = eig(layer2); 
          [EvL1,I1] = sort(diag(DL1));  VL1 = VL1(:,I1); 
          [EvL2,I2] = sort(diag(DL2));  VL2 = VL2(:,I2); 

          EigL1_trial_array(:,jtrial) = EvL1;
          EigL2_trial_array(:,jtrial) = EvL2;
          
          ind1 = 2*((inet1-1)*4 + inet2)-1
          f1=figure(ind1); 
          subplot(5,5,jdens)
          plot(EigL1_trial_array(:,jtrial)); 
          title(['idens = ',num2str(jdens)], 'FontWeight','normal')
          ylabel('Eigenval')
          xlabel('index')
          hold on
          axis square
          drawnow
          
          ind2 = 2*((inet1-1)*4 + inet2)
          f2 = figure(ind2); 
          subplot(5,5,jdens)
          plot(EigL2_trial_array(:,jtrial)); 
          title(['idens = ',num2str(jdens)], 'FontWeight','normal')
          ylabel('Eigenval')
          xlabel('index')
          hold on
          axis square
          drawnow
          
        end
    
    end
    
    % assign the density array
    EigL1_dens_arr(:,jdens) = mean(EigL1_trial_arr,2); 
    EigL2_dens_arr(:,jdens) = mean(EigL2_trial_arr,2); 
   
    
    
    if Adj_plot ==1 && fig_save ==1
    saveas(gcf,fullfile(fig_dir,['diff_SupA_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'idens=',num2str(jdens),'.png']))
    close(gcf)
    end
   
    if Eig_plot ==1 && fig_save ==1
    saveas(f1,fullfile(fig_dir,['EigL1_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.png']) )
    saveas(f2,fullfile(fig_dir,['EigL2_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.png']) )
    end
    
    
    
    
    end
    
    
end
end