close all; 
% dirname = 'Full_data_mat'; 
fig_opt =1; 
ntrials = 50; 
N_rot = 18;
sv = linspace(0,1,N_rot); 
fig_dir = 'Figures_Layer_Rotation'; 

inet1 = 1;
inet2 = 2;
itrial =1; 

for inet1=1:4
for inet2=1:4

dirname = ['trial=',num2str(itrial)]; 
% load mat file
load(fullfile(dirname,['LayerRotation_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'trial=',num2str(itrial),'_original.mat']))

mov_obj = VideoWriter(['network_evol_net1=',num2str(inet1),'net2=',num2str(inet2),'.avi']);

for js =1:size(sv,2)
%     figure;
    Duplex = load(fullfile(dirname,['Duplex_inet1=',num2str(inet1),...
                 '_inet2=',num2str(inet2),'_s=',num2str(js), '.csv']));
    
   imagesc(Duplex); 
   axis square
   xlabel('nodes')
   ylabel('nodes'); 
   set(gca,'fontsize',14)
   colorbar

   drawnow
   open(mov_obj) 
   frame3 = getframe(gcf);
   writeVideo(mov_obj, frame3)
        
    
end
close(mov_obj)

% scatter((alignment),E2,'.'); hold on; 
% box on
% axis square
% xlabel('alignment')
% ylabel('energy'); 
% set(gca,'fontsize',14)
% drawnow 
    end
end