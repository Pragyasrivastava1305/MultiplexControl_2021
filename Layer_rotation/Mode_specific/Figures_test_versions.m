close all; 
clear all;
% dirname = 'Full_data_mat'; 
% fig_opt = 1; 
% options of figures 
E_vs_align_scatter_opt = 0; 
nchannels_vs_align_opt = 0; 
E_vs_s_opt = 0;  
align_vs_s_opt = 0; 
align_vs_s_trial_opt = 0;
clen =16;
% colors 
% yellow1 = [204, 153, 0]/255;  yellow2 = [255, 204, 153]/255; 
% mauve1  = [204,51,153]/255;    mauve2 = [255,153,255]/255; 
% green1 = [0, 153,51]/255;     green2 = [0,255,205]/255;
% blue1 = [102,0,255]/255;     blue2 = [0,204,255]/255; 
% 
% col1 = [linspace(mauve1(1),mauve2(1),clen)', linspace(mauve1(2),mauve2(2),clen)', linspace(mauve1(3),mauve2(3),clen)'];
% col2 = [linspace(green1(1),green2(1),clen)', linspace(green1(2),green2(2),clen)', linspace(green1(3),green2(3),clen)'];
% col3 = [linspace(yellow1(1),yellow2(1),clen)', linspace(yellow1(2),yellow2(2),clen)', linspace(yellow1(3),yellow2(3),clen)'];
% col4 = [linspace(blue1(1),blue2(1),clen)', linspace(blue1(2),blue2(2),clen)', linspace(blue1(3),blue2(3),clen)'];
% 
% bcol1 = col4(3,:); 
% mcol1 = col1(2,:);
% ycol1 = col3(1,:); 
% gcol1 = col2(1,:);  


ntrials = 50; 
N_rot = 18;
sv = linspace(0,1,N_rot); 
fig_dir = 'Figures_Layer_Rotation'; 

E_cutoff = 0.01; 
E_tol = 10e-3;
inet1 =2;
inet2 =2;

for inet1 = 1:4
    

for inet2 = 1:4
    figure;
    
% assign colors based on the topology of second layer
        if inet2 ==1
%               col_sch = col1;
              col_one = [95 158 160]/255;
        elseif inet2 ==2
%               col_sch = col2;
              col_one = [186,181,147]/255;
        elseif inet2 == 3
%               col_sch = col3;
              col_one = [167,196,139]/255;
        else
%               col_sch = col4;
              col_one = [231,130,162]/255; 
        end
      

[inet1 inet2]
E_all_trials = zeros(N_rot,ntrials);
align_all_trials = zeros(N_rot,ntrials);
% itrial = 1;

for itrial =1:ntrials

% assign directory
dirname = ['trial=',num2str(itrial)]; 

% load mat file
load(fullfile(dirname,['LayerRotation_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'trial=',num2str(itrial),'_original.mat']))
sd_chk(itrial,:) = sd;

for js =1:size(sv,2) 

E_all_trials(:,itrial) = E2;
align_all_trials(:,itrial) = (alignment);
    
% reshape the 'projection of optim_u2 along the modes of first layer'
proj_reshape = reshape(proj_array(:,js),[size(optim_u2,1),size(optim_u2,2)]); 
    
% calculate the energy carried by each channel 
for ichannel =1:N
     E_channel(ichannel, js) = trapz(tarray,proj_reshape(ichannel,:).^2);
end

% get the maximum energy and index of corresponding mode
[max_E(js), channel_id(js)] = max(E_channel(:,js)); 

% find the number of activated channels
%     n_channel(itrial,js) = size( find( E_channel(:,js)>= 0.01*max_E(js) ),1); 
n_channel(js,itrial) = size( find( E_channel(:,js)> E_tol),1); 
    
end

if E_vs_align_scatter_opt ==1
plot(abs(align_all_trials(:,itrial)),E_all_trials(:,itrial),'-','color',col_one,'linewidth',1,'MarkerEdgeColor')
hold on
% splot = scatter(abs(alignment),E2,'o'); hold on; 
box on
axis square
xlabel('Alignment')
ylabel('Energy'); 
set(gca,'fontsize',14)
drawnow 
end


for js =1:size(sv,2)
    plot(align_all_trials(js,:),'-o'); hold on 
    set(gca,'fontsize',14)
    xlabel('trials')
    title(['inet1=',num2str(inet1),'inet2=',num2str(inet2)])
    drawnow
end




% if nchannels_vs_align_opt ==1
% % figure(inet1)
% plot(abs(align_all_trials(:,itrial)), n_channel(:,itrial),'-','color',col_one,'linewidth',1,'MarkerEdgeColor','k')
% hold on
% xlabel('Alignment')
% ylabel('Number of channels'); 
% % ylim([0 250])
% set(gca,'fontsize',14)
% % title(['inet1=',num2str(inet1)])
% title(['inet1=',num2str(inet1),' inet2=',num2str(inet2)], 'FontWeight','normal')
% % title(['inet1=',num2str(inet1),'inet2=',num2str(inet2)])
% axis square
% drawnow
% end


if align_vs_s_trial_opt ==1
% figure(inet1)
 plot(sv', align_all_trials(:,itrial),'sq','color',col_one,'linewidth',1.5,'MarkerEdgeColor',[0.7 0.7 0.7])

hold on
xlabel('s')
ylabel('Alignment'); 
% ylim([0 250])
set(gca,'fontsize',14)
% title(['inet1=',num2str(inet1)])
title(['inet1=',num2str(inet1),' inet2=',num2str(inet2)], 'FontWeight','normal')
axis square
drawnow
end




end


if nchannels_vs_align_opt ==1
plot(mean(abs(align_all_trials),2), mean(n_channel,2),'-s','color',col_one,'MarkerSize',10,'LineWidth',2,'MarkerEdgeColor',col_one)
hold on
axis square
ylim([1 100])
grid on
xlabel('Alignment')
ylabel('Number of excited channels')
set(gca,'FontSize',14)
legend('ER','WS','BA','RG')
drawnow
end


if E_vs_align_scatter_opt ==1
plot(mean(abs(align_all_trials),2), mean(E_all_trials,2),'-o','color','col_one','linewidth',2,'MarkerEdgeColor','col_one','MarkerSize',10)
hold on
axis square
ylim([1 100])
grid on
xlabel('Alignment')
ylabel('Number of excited channels')
set(gca,'FontSize',14)
legend('ER','WS','BA','RG')
drawnow
end



if align_vs_s_opt ==1
subplot(1,3,1)
err_align = std(align_all_trials,0,2);
% plot(s,abs(alignment),'LineWidth',1.5); hold on; 
errorbar(sv', mean(align_all_trials,2),err_align,'-s','MarkerSize',10,'LineWidth',1.5)
hold on
xlabel('s')
ylabel('alignment'); 
set(gca,'fontsize',14)
% title(['inet1=',num2str(inet1),'inet2=',num2str(inet2)])
title(['inet1=',num2str(inet1)])
axis square
drawnow
end

if E_vs_s_opt ==1
subplot(1,3,2)
err_E = std(E_all_trials,0,2);
errorbar(sv', mean(E_all_trials,2),err_E,'-s','MarkerSize',10,'LineWidth',1.5)
hold on
xlabel('s')
ylabel('energy'); 
ylim([0 250])
set(gca,'fontsize',14)
title(['inet1=',num2str(inet1)])
% title(['inet1=',num2str(inet1),'inet2=',num2str(inet2)])
axis square
drawnow
end


% saveas(gcf, fullfile(fig_dir,['Alignment_vs_s_all_trials_inet1=',num2str(inet1),'inet2=',num2str(inet2), '.png'])) ; 
% saveas(gcf, fullfile(fig_dir,['nchannels_vs_Align_all_trials_inet1=',num2str(inet1),'inet2=',num2str(inet2), '.png'])) ;
end
        

 
end
