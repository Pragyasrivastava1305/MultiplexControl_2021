close all; 
clear all;

ttrials = 100; 
N_rot = 18;
sv = linspace(0,1,N_rot); 
fig_dir = 'Figures_Layer_Rotation'; 

inet1 =3;
inet2 =3;

mk_array = {'o','sq','<','d'};

for inet1 =1:4
for inet2 =1:4
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
E_all_trials = zeros(N_rot,ttrials);
align_all_trials = zeros(N_rot,ttrials);

for jtrial =1:ttrials

% assign directory
dirname = ['trial=',num2str(jtrial)]; 

% load mat file
load(fullfile(dirname,['LayerRotation_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'trial=',num2str(jtrial),'_original.mat']))

for js =1:size(sv,2) 

align_all_trials(:,jtrial) = (alignment);
E_all_trials(:,jtrial) = E2;
    
end
end

[xfull, Ifull] = sort((align_all_trials(:)));
yt = E_all_trials(:);
yfull = yt(Ifull); 

% Ilow = find(abs(yfull) );
% xl = round(mean(abs(xfull(Ilow))),2);

plot(xfull,yfull,'s','color',col_one,'MarkerFaceColor',col_one,'MarkerEdgeColor',col_one ); hold on 

[xvec,yvec,min_vec,max_vec] = binned_vectors((align_all_trials), E_all_trials, 29,[-1 1]);

col_one ='k'; 
plot(xvec,yvec, '-', 'Marker',mk, 'Color','k', 'LineWidth',2, 'MarkerSize',10,...
           'MarkerFaceColor',col_one,'MarkerEdgeColor',col_one); hold on

% jbfill(xvec, movmean(min_vec,3), movmean(max_vec,3),col_one,col_one, 1, 0.15)
       

hold on; 
set(gca,'fontsize',14)
xlabel('alignment')
ylabel('Energy')
grid on
drawnow

saveas(gcf, fullfile(fig_dir,['E_vs_Align_scatter_inet1=',num2str(inet1),'inet2=',num2str(inet2), '.pdf'])) 
end
% saveas(gcf, fullfile(fig_dir,['E_vs_Align_scatter_inet1=',num2str(inet1), '.pdf'])) 

% legend('ER', 'WS', 'BA', 'RG')
end
