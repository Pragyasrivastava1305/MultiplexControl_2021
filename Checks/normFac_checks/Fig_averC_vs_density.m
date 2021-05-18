clear all
close all;

N = 100;
num_trial = 2;
ndens = 11; 
ntop = 4;

% define colors
def_colors;
inet2 = 4;
for inet1 =1:ntop
% inet2 =1
%  
for inet2  = 1:ntop
    % assign colors according to the topology of first layer
    if inet1 ==1
        col_sch = col1;
        col_one = bcol1;
    elseif inet1 ==2
        col_sch = col2;
        col_one = ycol1;
    elseif inet1 == 3
        col_sch = col3;
        col_one = gcol1;
    else
        col_sch = col4;
        col_one = mcol1; 
    end
          
    if inet2 ==1
        mk = 'o';
    elseif inet2 ==2
        mk = 'sq';
    elseif inet2 == 3
        mk = 'd';    
    else
        mk = 'o'; 
    end
          
    mk_array = {'o','sq','<','d'};
        
    E1_aver = zeros(N,ndens); 
    E2_aver = zeros(N,ndens);
    rho2_array = zeros(ndens,num_trial); 
    
    for jtrial =1:num_trial
        dirname = ['norm_check/trial=',num2str(jtrial),'/inet1=',num2str(inet1),'_inet2=',num2str(inet2)]; 

        % load mat file
        load(fullfile(dirname,['topodensity_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'itrial=',num2str(jtrial),'.mat']))


        E1_aver = E1_aver + E1_array; 
        E2_aver = E2_aver + E2_array; 
        rho_vec(jtrial,inet1) = rho1;
        rho2_array(:,jtrial) =  rho2; 
    end
    
    % average over trials
    E1_aver = E1_aver/num_trial; 
    E2_aver = E2_aver/num_trial;  
    xvec = mean(rho2_array,2)';
    
%    % define average control energy and maximum control energy: will have
%    % same length as density vec
     averC_L1 = mean(E1_aver); 
     averC_L2 = mean(E2_aver); 
     MaxC_L1 = max(E1_aver); 
     MaxC_L2 = max(E2_aver); 
     
     Evec1 = averC_L2; 
     Evec2 = MaxC_L2; 
     
     
     figure(1)
     subplot(2,4, inet2)
     plot(xvec, Evec1,'o-','color',col_one,'MarkerFaceColor',col_one, 'MarkerEdgeColor',[1 1 1],'linewidth',2,...
                               'MarkerSize',8)
                          
     hold on; 
     set(gca,'fontsize',14,'fontweight','bold')
     % xlabel('\Delta \rho')  
     % ylim([0.82 0.88])
     xlim([0 1])
     axis square
     drawnow
        
     subplot(2,4, inet2 + ntop)
     plot(xvec, Evec2,'o-','color',col_one,'MarkerFaceColor',col_one, 'MarkerEdgeColor',[1 1 1],'linewidth',2,...
                              'MarkerSize',8)
                          
     hold on; 
     set(gca,'fontsize',14,'fontweight','bold')
     xlim([0 1])
     % ylim([ 60 125])
     % xlabel('\Delta \rho')
     %ylabel('Average control energy for layer 2','FontWeight', 'normal')
     axis square
% %    leg  = legend('ER','WS','BA','RG'); 
%      title(leg,'layer 1')
     drawnow
     
     figure(2)
      plot(max_eig_G1); hold on;  drawnow
                 
end 
       
end 
% % 
% % 
% % saveas(gcf,fullfile('figures',['L2_aver&mod_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.fig']))
% % saveas(gcf,fullfile('figures',['L2_aver&mod_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.png']))

