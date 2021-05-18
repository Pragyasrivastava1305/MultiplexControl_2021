clear all
close all;
ind=25; 

% define colors
ntop = 4;
strt_ind = 5;
dirname = 'Full_data_mat_rho=05'; 

for jnet1 =1:ntop
for jnet2  = 1:ntop

    cols_n_markers;
        
    load(fullfile(dirname,['EigenEnergy_inet1=',num2str(jnet1),'_inet2=',num2str(jnet2),'_original.mat']))

    rho_vec(jnet1,jnet2) = rho1;
        
        % initiate arrays to store averages over samplees 
        E1_aver = zeros(N,size(loop_vec,2)); 
        E2_aver = zeros(N,size(loop_vec,2)); 
        xi_aver = zeros(N,size(loop_vec,2)); 
        mu_aver = zeros(N,size(loop_vec,2)); 
        
        for idens = 1:size(loop_vec,2)
            
            E1_mat =  reshape(E1_array(:,idens),[N, ntrial]); 
            E2_mat =  reshape(E2_array(:,idens),[N, ntrial]); 
            xi_mat =  reshape(xi_array(:,idens),[N, ntrial]); 
            mu_mat =  reshape(mu_array(:,idens),[N, ntrial]); 
            
            E1_aver(:,idens) = mean(E1_mat,2)';
            E2_aver(:,idens) = mean(E2_mat,2)';
            xi_aver(:,idens) = mean(xi_mat,2)';
            mu_aver(:,idens) = mean(mu_mat,2)';
            dspec(idens) = sqrt(sum(( xi_aver(:,idens) - mu_aver(:,idens)).^2)/N ) ; 
           
            % calculate average control energy and min/max 
            av_E_L1(idens,jnet1) = mean(E1_aver(:,idens)); 
            [min_E_L1(idens,jnet1),imax_L1(idens,jnet1)] = min(E1_aver(:,idens)); 
           [max_E_L1(idens,jnet1),imin_L1(idens,jnet1)] = max(E1_aver(:,idens)); 
            
            av_E_L2(idens,jnet1) = mean(E2_aver(:,idens)); 
           [min_E_L2(idens,jnet1),imax_L2(idens,jnet1)] = min(E2_aver(:,idens)); 
           [max_E_L2(idens,jnet1),imin_L2(idens,jnet1)] = max(E2_aver(:,idens)); 
            
        end  
       
        
        Evec1 = av_E_L2; 
        Evec2 = max_E_L2; 
        
        figure(1)
        subplot(2,4, jnet2)
        plot(mean(rho2(1:end,:),2), Evec1(1:end,jnet1),'Marker',mk,'color',col_one,...
                          'MarkerFaceColor',col_one, 'MarkerEdgeColor',[1 1 1],'linewidth',2,...
                              'MarkerSize',8)
                          
        hold on; 
        set(gca,'fontsize',14,'fontweight','bold')
        %xlabel('\Delta \rho')  
        %ylim([0.83 0.88])
        ylim([200 350])
        xlim([0 1])
        axis square
        grid on
        drawnow
        
        subplot(2,4, jnet2 + ntop)
        plot(mean(rho2(1:end,:),2), Evec2(1:end,jnet1),'Marker', mk,'color',col_one...
                      ,'MarkerFaceColor',col_one, 'MarkerEdgeColor',[1 1 1],'linewidth',2,...
                              'MarkerSize',8)
                          
        hold on; 
        set(gca,'fontsize',14,'fontweight','bold')
        xlim([0 1])
        ylim([100 5000])
        set(gca,'yscale','log')
        axis square
        grid on
%       leg  = legend('ER','WS','BA','RG'); 
%       title(leg,'layer 1')
        drawnow
        
        
end 
      
end 

% saveas(gcf,fullfile('figures',['L2_aver&mod_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.fig']))
% saveas(gcf,fullfile('figures',['L2_aver&mod_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.png']))
%          









