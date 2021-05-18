clear all
close all;

N = 100;
num_trial = 40;
ndens = 25; 
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
          
%     if inet2 ==1
%         mk = 'o';
%     elseif inet2 ==2
%         mk = 'sq';
%     elseif inet2 == 3
%         mk = 'd';    
%     else
%         mk = 'o'; 
%     end
    
    if inet2 ==1
              mk = 'o';
        elseif inet2 ==2
              mk = 'sq';
        elseif inet2 == 3
              mk = '<';    
        else
              mk = 'd'; 
    end
    mk_array = {'o','sq','<','d'};
        
    E1_aver = zeros(40,ndens); 
    E2_aver = zeros(40,ndens);
    E1_max = zeros(40, ndens); 
    E2_max = zeros(40, ndens); 
    rho2_array = zeros(ndens,num_trial); 
    
    for jtrial =1:40
        dirname = ['trial=',num2str(jtrial),'/inet1=',num2str(inet1),'_inet2=',num2str(inet2)]; 

        % load mat file
        load(fullfile(dirname,['topodensity_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'itrial=',num2str(jtrial),'.mat']))


        E1_aver(jtrial,:) = mean(E1_array); 
        E2_aver(jtrial,:) = mean(E2_array); 
        E1_max(jtrial,:) = max(E1_array); 
        E2_max(jtrial,:) = max(E2_array); 
        
        rho_vec(jtrial,inet1) = rho1;
        rho2_array(:,jtrial) =  rho2; 
    end
      
    xvec = mean(rho2_array,2)';
    
    E1_quant = E2_aver; 
    E2_quant = E2_max; 
    
%    % define average control energy and maximum control energy: will have
%    % same length as density vec
    mean_E1 = mean(E1_quant); 
    std_E1 = std(E1_quant); 
    up_c1 = mean_E1 + 0.5*std_E1; 
    lo_c1 = mean_E1 - 0.5*std_E1; 
    
    mean_E2 = mean(E2_quant); 
    std_E2 = std(E2_quant); 
    up_c2 = mean_E2 + 0.5*std_E2; 
    lo_c2 = mean_E2 - 0.5*std_E2;  
    
     figure(1)
     subplot(2,4, inet2)
     plot(xvec, mean_E1,'-','marker',mk,'color',col_one,'MarkerFaceColor',col_one, 'MarkerEdgeColor',[1 1 1],'linewidth',1.5,...
                               'MarkerSize',8)
%      hold on
%      errorbar(xvec, mean_E1, std_E1,'linestyle','none','color',[0.8 0.8 0.8],'linewidth',1.5)
                               
    %shade(xvec,lo_c1,xvec,up_c1,'Filltype',[1 2])
                          
     hold on; 
     set(gca,'fontsize',14,'fontweight','bold')
     % xlabel('\Delta \rho')  
     ylim([60 140])
     xlim([0 1])
     axis square
     drawnow
        
     subplot(2,4, inet2 + ntop)
     plot(xvec, mean_E2,'-','marker',mk,'color',col_one,'MarkerFaceColor',col_one, 'MarkerEdgeColor',[1 1 1],'linewidth',1.5,...
                               'MarkerSize',8)
     %subplot(1,2,2)
     
     hold on; 
%      errorbar(xvec, mean_E2,std_E2,'linestyle','none','color',[0.8 0.8 0.8],'linewidth',1.5)
     
     set(gca,'fontsize',14,'fontweight','bold')
     xlim([0 1])
     % ylim([ 60 125])
     % xlabel('\Delta \rho')
     %ylabel('Average control energy for layer 2','FontWeight', 'normal')
     axis square
% %    leg  = legend('ER','WS','BA','RG'); 
%      title(leg,'layer 1')
     drawnow
                 
end 
       
end 
% % 
% % 
% % saveas(gcf,fullfile('figures',['L2_aver&mod_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.fig']))
% % saveas(gcf,fullfile('figures',['L2_aver&mod_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.png']))

