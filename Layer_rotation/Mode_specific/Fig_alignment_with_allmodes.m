% This code checks that the rotation of the target layer modes relative to
% the first layer induces alignment changes relative to only those axes
% that form the rotation plane. alignment with other modes remain
% unchanged. 


close all; 
clear all;
% dirname = 'Full_data_mat'; 
% fig_opt = 1; 
% options of figures 

num_trials = 100; 
N_rot = 18;
sv = linspace(0,1,N_rot); 
fig_dir = 'Figures_Layer_Rotation/check_alignment_with_all_modes'; 

yellow2 = [255, 204, 153]/255;
blue2 = [0,204,255]/255; 

col1 = [linspace(yellow2(1),blue2(1),num_trials)',...
        linspace(yellow2(2),blue2(2),num_trials)',...
        linspace(yellow2(3),blue2(3),num_trials)'];


for inet1 = 1:4
    for inet2 = 1:4
    figure;
    
    align_all_trials = zeros(N_rot,num_trials);
    % itrial = 1;

    for jtrial =1:num_trials
        % assign directory
        if mod(jtrial,25)==0
        dirname = ['trial=',num2str(jtrial)]; 

        % load mat file
        load(fullfile(dirname,['LayerRotation_inet1=',num2str(inet1),...
                '_inet2=',num2str(inet2),'trial=',num2str(jtrial),'_original.mat']))

        subplot(2,2,jtrial/25)
        for js =1:size(sv,2) 
%             plot(1:N, abs(alignment),'-','color',[col1(jtrial,1),col1(jtrial,2),col1(jtrial,3)], 'MarkerEdgeColor','k','LineWidth',1); hold on 
%             plot(1:N, abs(alignment),'-', 'MarkerfaceColor','w','LineWidth',1); hold on 
            imagesc(1:N,sv,abs(alignment)')
            set(gca,'fontsize',14)
            xlabel('Eigenmode #')
            ylabel('alignment')
            title(['trial=',num2str(jtrial)],'FontWeight','normal')
            drawnow
        end
        end
        saveas(gcf, fullfile(fig_dir,['nchannels_vs_Align_scatter_inet1=',num2str(inet1),'inet2=',num2str(inet2), '.png'])) 
        
        
    end
    
    end
    
end
