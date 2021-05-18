% code to indentify the input eigenmode 
close all;  clear all;

t_trials = 100; 
N_rot = 18;
sv = linspace(0,1,N_rot); 
fig_dir = 'Figures_Layer_Rotation'; 
E_tol = 10e-3;

col1 = [134, 207, 190]/255; 
col2 = [ 251, 182, 209]/255; 
col3 = [0.4 0.4 0.4]; 

v = VideoWriter('test.avi'); 

for inet1 = 1:4
for inet2 = 1:4
    
% figure(inet1)
% call code to get colors and markers
cols_n_markers; 
[inet1 inet2]

align_all_trials = zeros(N_rot,t_trials);
jtrial = 1; 

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

        % get the maximum energy and index of corresponding mode
        [max_E(jtrial, js), channel_id(jtrial, js)] = max(E_channel(:,js)); 
        
        % calculate the fraction of energy carried by first two modes 
        Efrac1(jtrial,js) = E_channel(100,js)/sum(E_channel(:,js)); 
        Efrac2(jtrial,js) = E_channel(99,js)/sum(E_channel(:,js)); 
        
        % calculate the fraction of energy carried by the maximum energy
        % mode
        Efrac0(jtrial,js) = max_E(jtrial,js)/sum(E_channel(:,js));
        

        % find the number of activated channels
        % n_channel(itrial,js) = size( find( E_channel(:,js)>= 0.01*max_E(js) ),1); 
        % n_channel(js,itrial) = size( find( E_channel(:,js)> E_tol),1); 
    
    end
     figure(1)
     subplot(4,4,(inet1-1)*4 + inet2)
     
     plot(abs(alignment(end,:)), Efrac1(jtrial,:),'o','MarkerFaceColor'...
                                            ,col1, 'MarkerEdgeColor', 'k');  
     box on
     xlim([0 1])
     ylim([0 1])
     set(gca,'fontsize',14)
     axis square
     drawnow
     
    hold on; 
    if inet1 ~=4 
        xticklabels([])
    end

    if inet2 ~=1
        yticklabels([])
    end
     
%      plot(abs(alignment(end,:)), Efrac2(jtrial,:),'o','MarkerFaceColor'...
%                                             ,col2, 'MarkerEdgeColor','k'); 

     
%      figure(2)
%      subplot(4,4,(inet1-1)*4 + inet2)
%      plot(abs(alignment(end-1,:)), Efrac2(jtrial,:),'o','MarkerFaceColor'...
%                                              ,col2, 'MarkerEdgeColor', 'k'); 
%       
%      box on
%      xlim([0 1])
%      ylim([0 1])
%      set(gca,'fontsize',14)
%      axis square
%      drawnow
%      
%      hold on; 
%      
%      if inet1 ~=4 
%         xticklabels([])
%      end
% 
%      if inet2 ~=1
%         yticklabels([])
%      end

     figure(3)
     subplot(4,4,(inet1-1)*4 + inet2)
     plot(abs(alignment(end,:)), Efrac1(jtrial,:) + Efrac2(jtrial,:),'o','MarkerFaceColor'...
                                            ,[0.8 0.8 0.8], 'MarkerEdgeColor', 'k','MarkerSize',6);  
     box on
     xlim([0 1])
     ylim([0 1])
     set(gca,'fontsize',14)
     axis square
     drawnow
     
     hold on; 
      if inet1 ~=4 
        xticklabels([])
    end

    if inet2 ~=1
        yticklabels([])
    end

end

end
end

