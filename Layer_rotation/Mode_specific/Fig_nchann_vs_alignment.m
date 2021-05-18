close all;  clear all;

t_trials = 100; 
N_rot = 18;
sv = linspace(0,1,N_rot); 
fig_dir = 'Figures_Layer_Rotation'; 
E_tol = 10e-3;

for inet1 =1:4
for inet2 =1:4
    
    % figure(inet1)
    % call code to get colors and markers
    cols_n_markers; 
        
    [inet1 inet2]
    align_all_trials = zeros(N_rot,t_trials);

    for jtrial =1:t_trials
        % assign directory
        dirname = ['trial=',num2str(jtrial)]; 

        % load mat file
        load(fullfile(dirname,['LayerRotation_inet1=',num2str(inet1),...
                '_inet2=',num2str(inet2),'trial=',num2str(jtrial),'_original.mat']))

        for js =1:size(sv,2) 

        align_all_trials(:,jtrial) = (alignment(end,:));
    
        % reshape the 'projection of optim_u2 along the modes of first layer'
        proj_reshape = reshape(proj_array(:,js),[size(optim_u2,1),size(optim_u2,2)]); 
    
        % calculate the energy carried by each channel 
        for ichannel =1:N
            E_channel(ichannel, js) = trapz(tarray,proj_reshape(ichannel,:).^2);
        end

        % find the number of activated channels
%       n_channel(itrial,js) = size( find( E_channel(:,js)>= 0.01*max_E(js) ),1); 
        n_channel(js,jtrial) = size( find( E_channel(:,js)> E_tol),1); 
    
        end
end

[xfull, Ifull] = sort(abs(align_all_trials(:)));
yt = n_channel(:);
yfull = yt(Ifull); 

% plot(xfull,yfull,'s','color',col_one,'MarkerFaceColor',col_one,'MarkerEdgeColor',col_one ); hold on 

% find the end-points of the x-axis
Ilow = find(abs(yfull) ==1);
xl = round(mean(abs(xfull(Ilow))),2);

[xvec,yvec,min_vec,max_vec] = binned_vectors(abs(align_all_trials), n_channel, 14,[0 xl]); 
% 
% col_one = 'k'; 
plot(xvec,yvec, '-', 'Marker',mk, 'Color',col_one, 'LineWidth',2, 'MarkerSize',14,...
            'MarkerFaceColor',col_one,'MarkerEdgeColor','w'); hold on

%plot(xvec,yvec, '-', 'Marker',mk, 'Color',col_one, 'LineWidth',2, 'MarkerSize',10,...
%            'MarkerFaceColor','k','MarkerEdgeColor','w'); hold on
%jbfill(xvec, movmean(min_vec,3), movmean(max_vec,3),col_one,col_one, 1, 0.15)
        

hold on; 
set(gca,'fontsize',14)
xlabel('alignment')
ylabel('Number of channels')
grid on
drawnow

% saveas(gcf, fullfile(fig_dir,['nchannels_vs_Align_scatter_inet1=',num2str(inet1),'inet2=',num2str(inet2), '.pdf'])) 
end
%saveas(gcf, fullfile(fig_dir,['nchannels_vs_Align_scatter_inet1=',num2str(inet1), '.pdf'])) 
%saveas(gcf, fullfile(fig_dir,['nchannels_vs_Align_all_trials_inet1=',num2str(inet1), '.pdf'])) 
% 
% 

end
legend('ER-ER', 'ER-WS', 'ER-BA', 'ER-RG', 'WS-ER', 'WS-WS','WS-BA', 'WS-RG', ...
            'BA-ER', 'BA-WS','BA-BA', 'BA-RG','RG-ER', 'RG-WS','RG-BA', 'RG-RG' )