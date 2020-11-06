close all; 
clear all;

t_trials = 100; 
N_rot = 18;
sv = linspace(0,1,N_rot); 
fig_dir = 'Figures_Layer_Rotation'; 
E_tol = 10e-3;
inet1 =4; 
inet2 =4;

for inet1 =1:4
for inet2 =1:4
figure(inet1)
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
      
        if inet1 ==1
              mk = 'o';
        elseif inet1 ==2
              mk = 'sq';
        elseif inet1 == 3
              mk = '+';    
        else
              mk = '<'; 
        end
        
        
[inet1 inet2]
align_all_trials = zeros(N_rot,t_trials);

for itrial =1:t_trials

% assign directory
dirname = ['trial=',num2str(itrial)]; 

% load mat file
load(fullfile(dirname,['LayerRotation_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'trial=',num2str(itrial),'_original.mat']))

for js =1:size(sv,2) 

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
plot(xvec,yvec, '-', 'Marker',mk, 'Color',col_one, 'LineWidth',2, 'MarkerSize',10,...
            'MarkerFaceColor',col_one,'MarkerEdgeColor',col_one); hold on
jbfill(xvec, movmean(min_vec,3), movmean(max_vec,3),col_one,col_one, 1, 0.15)
        

hold on; 
set(gca,'fontsize',14)
xlabel('alignment')
ylabel('Number of channels')
grid on
drawnow

% saveas(gcf, fullfile(fig_dir,['nchannels_vs_Align_scatter_inet1=',num2str(inet1),'inet2=',num2str(inet2), '.pdf'])) 
end
saveas(gcf, fullfile(fig_dir,['nchannels_vs_Align_scatter_inet1=',num2str(inet1), '.pdf'])) 
% % saveas(gcf, fullfile(fig_dir,['nchannels_vs_Align_all_trials_inet1=',num2str(inet1), '.pdf'])) 
% 
% 
% % legend('ER', 'WS', 'BA', 'RG')
end
