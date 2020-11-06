% code to check if 
close all; 
clear all;
% dirname = 'Full_data_mat'; 
% fig_opt = 1; 
% options of figures 

t_trials = 100; 
N_rot = 18;
sv = linspace(0,1,N_rot); 
fig_dir = 'Figures_Layer_Rotation'; 

E_cutoff = 0.01; 
E_tol = 10e-3;


for inet1 = 1:4
    for inet2 = 1:4
    figure;
    
    align_all_trials = zeros(N_rot,t_trials);
    % itrial = 1;

    for jtrial =1:t_trials
        % assign directory
        dirname = ['trial=',num2str(jtrial)]; 

        % load mat file
        load(fullfile(dirname,['LayerRotation_inet1=',num2str(inet1),...
                '_inet2=',num2str(inet2),'trial=',num2str(jtrial),'_original.mat']))

        for js =1:size(sv,2) 
            align_all_trials(:,jtrial) = abs(alignment(end,:));
        end
        
        plot(sv',align_all_trials(:,jtrial),'-s','color',[0.5 0.5 0.5], 'MarkerEdgeColor','k','LineWidth',1); hold on 
        set(gca,'fontsize',14)
        xlabel('s')
        ylabel('alignment')
        title(['inet1=',num2str(inet1),' inet2=',num2str(inet2)],'FontWeight','normal')
        drawnow
        
    end

%     
%     for js =1:size(sv,2)
%         plot(align_all_trials(js,:),'-o'); hold on 
%         set(gca,'fontsize',14)
%         xlabel('trials')
%         title(['inet1=',num2str(inet1),'inet2=',num2str(inet2)])
%         drawnow
%     end

%     saveas(gcf, fullfile(fig_dir,['align_check_inet1=',num2str(inet1),'inet2=',num2str(inet2), '.png'])) ; 
   % saveas(gcf, fullfile(fig_dir,['align_vs_s_inet1=',num2str(inet1),'inet2=',num2str(inet2), '.png'])) ; 
    end
   
    
    
end
