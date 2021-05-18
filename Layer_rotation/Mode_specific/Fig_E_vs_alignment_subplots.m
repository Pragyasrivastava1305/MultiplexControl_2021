% code to plot E vs alignment, with curves for each plotted in a different
% subplot

close all; 
%clear all;

ttrials = 100; 
N_rot = 18;
sv = linspace(0,1,N_rot); 
fig_dir = 'Figures_Layer_Rotation'; 

for inet1 =1:4
for inet2 =1:4

    figure(1)
    % call this code to get colors and markers
    cols_n_markers; 
    [inet1 inet2]
    
    E_all_trials = zeros(N_rot,ttrials);
    align_all_trials = zeros(N_rot,ttrials);

    for jtrial =1:ttrials
        % assign directory
        dirname = ['trial=',num2str(jtrial)]; 
    
        % load mat file
        load(fullfile(dirname,['LayerRotation_inet1=',num2str(inet1)...
            ,'_inet2=',num2str(inet2),'trial=',num2str(jtrial),'_original.mat']))

        for js =1:size(sv,2) 
            align_all_trials(:,jtrial) = (alignment(end,:));
            E_all_trials(:,jtrial) = E2;
        end
    end

    subplot(4,4,(inet1-1)*4 + inet2)
    [xfull, Ifull] = sort(abs(align_all_trials(:)));
    yt = E_all_trials(:);
    yfull = yt(Ifull); 

    plot(xfull,yfull,'s','color',col_one,'MarkerFaceColor',col_one,'MarkerEdgeColor',col_one ); hold on 
    %plot(xfull, yfull)

    xl = min((xfull)); 
    xh = max((xfull)); 
    [xvec,yvec,min_vec,max_vec] = binned_vectors((xfull), yfull, 29,[xl xh]);

    col_one ='k'; 
    plot(xvec,movmean(yvec,2),'LineWidth',2, 'Color', col_one); hold on

    ylim([0 300])

    hold on; 
    set(gca,'fontsize',14)
    set(gca,'yscale','log')

    if inet1 ~=4 
        xticklabels([])
    end

    if inet2 ~=1
        yticklabels([])
    end

    grid on
    drawnow

    %saveas(gcf, fullfile(fig_dir,['E_vs_Align_scatter_inet1=',num2str(inet1),'inet2=',num2str(inet2), '.pdf'])) 
end
% saveas(gcf, fullfile(fig_dir,['E_vs_Align_scatter_inet1=',num2str(inet1), '.pdf'])) 

% legend('ER', 'WS', 'BA', 'RG')
end
