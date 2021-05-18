% this code produces plots of Ei's vs eigenvalues plots, where 
% Ei's are averaged over all the trials. For simplicity, all the densities
% are plotted with one color. 

clear all
close all;
N =100; 
ndens = 25;
ntop =4;
num_trial =40;
ysc = 2;  % if ysc =1, linear scale, if ysc =2, log scale.

% define colors and markers
mk_arr = ["o", "sq", "<", "d"]; 
c_arr = ["r","b","g","k"]; 
def_colors;

% define plot options
fig_dir = 'Figures/E_vs_evals';
opt_trial_plts = 0; 
opt_aver_plts = 1;
plt_save = 0;


for inet1 =1:ntop
% % %inet2 =1
% 
for inet2  = 1:ntop
    
    % assign colors according to the topology of first layer
    if inet2 ==1
        col_one = bcol1;
    elseif inet2 ==2
        col_one = ycol1;
    elseif inet2 == 3
        col_one = gcol1;
    else
        col_one = mcol1; 
    end
    
    
    E1_aver = zeros(N,ndens); 
    E2_aver = zeros(N,ndens);
    EigL1_aver = zeros(N,ndens); 
    EigL2_aver = zeros(N,ndens); 
    rho2_array = zeros(ndens,num_trial);
    jtrial =1;  
    
    for jtrial =1: num_trial
        dirname = ['trial=',num2str(jtrial),'/inet1=',num2str(inet1),'_inet2=',num2str(inet2)]; 
    
        % load mat file
        load(fullfile(dirname,['topodensity_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'itrial=',num2str(jtrial),'.mat']))
    
        E1_aver = E1_aver + E1_array; 
        E2_aver = E2_aver + E2_array; 
        rho_vec(jtrial,inet1) = rho1;
        rho2_array(:,jtrial) =  rho2; 
  
        EigL1_dens_arr = zeros(N, ndens);
        EigL2_dens_arr = zeros(N, ndens);
        
        for jdens = 1: ndens
            SupraA = load(fullfile(dirname,['SupA_inet1=',num2str(inet1),'_inet2=',...
                                                 num2str(inet2),'_trial=',num2str(jtrial),'_dens=',num2str(jdens),'.csv']));  
            layer1 = SupraA(1:N,1:N); 
            layer2 = SupraA(N+1:2*N, N+1:2*N); 

            [VL1,DL1] = eig(layer1); 
            [VL2,DL2] = eig(layer2); 
            [EvL1,I1] = sort(diag(DL1));  VL1 = VL1(:,I1); 
            [EvL2,I2] = sort(diag(DL2));  VL2 = VL2(:,I2); 

            EigL1_dens_arr(:,jdens) = EvL1;
            EigL2_dens_arr(:,jdens) = EvL2;
        end
        EigL1_aver = EigL1_aver + EigL1_dens_arr; 
        EigL2_aver = EigL2_aver + EigL2_dens_arr; 
    end
    
    % average over trials
    E1_aver = E1_aver/num_trial; 
    E2_aver = E2_aver/num_trial;  
    EigL1_aver = EigL1_aver/num_trial; 
    EigL2_aver = EigL2_aver/num_trial; 
    
    if opt_aver_plts ==1
        for jdens = 1: ndens
            if mod(jdens,5) ==0
            f1 = figure(1); 
            subplot(1,4,inet1)
           
                        
            plot(EigL1_aver(:,jdens), E1_aver(:,jdens),'o',...
                        'MarkerFaceColor',col_one,'MarkerEdgeColor',col_one); hold on; drawnow
            
            set(gca,'yscale','linear');  set(gca,'Fontsize',14,'fontweight','bold')
            xlabel('\xi_\alpha'); ylabel('E_\alpha')
            grid on; axis square;  box on

            
            f2 = figure(2); 
            subplot(1,4,inet1)
       
            plot(EigL2_aver(:,jdens), E2_aver(:,jdens),'o',...
                        'MarkerFaceColor',col_one,'MarkerEdgeColor',col_one); hold on; drawnow
            
            set(gca,'yscale','linear'); 
            set(gca,'Fontsize',14,'fontweight','bold')
            xlabel('\mu_\beta'); ylabel('E_\beta')
            grid on ; axis square; box on
            
            end
        end
        
        if plt_save ==1
           saveas(f1,fullfile(fig_dir,['L1_E_vs_eval_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.fig']))
           saveas(f1,fullfile(fig_dir,['L1_E_vs_eval_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.png']))
           saveas(f2,fullfile(fig_dir,['L2_E_vs_eval_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.fig']))
           saveas(f2,fullfile(fig_dir,['L2_E_vs_eval_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.png']))
        end
    end
    
end 
      
end 
% % 
% % 
% % saveas(gcf,fullfile('figures',['L2_aver&mod_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.fig']))
% % saveas(gcf,fullfile('figures',['L2_aver&mod_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.png']))

