close all; 
dirname = 'trial=5/inet1=1_inet2=2'; 
fig_opt =1; 
itrial =5; 
inet1 = 1; 
inet2 = 2;

% load mat file
SupA = load(fullfile(dirname,['SupA_inet1=',num2str(inet1),'_inet2=',...
                                         num2str(inet2),'_trial=',num2str(itrial),'_dens=',num2str(idens),...
                                                                                '.csv'])); 

fig_dir = 'Figs_projection_analysis'; 

idens = 25; 

%% 
% load duplex adjacency 
SupA = load(fullfile(dirname,['SupA_inet1=',num2str(inet1),'_inet2=',...
                                         num2str(inet2),'_trial=',num2str(itrial),'_dens=',num2str(idens),...
                                                                                '.csv'])); 

% imagesc(SupA) 
% extract layers 
first_layer = SupA(1:N,1:N); 
second_layer =  SupA(N+1:2*N, N+1:2*N); 

% eigen-decomposition of the two layers
[VL1,EL1] = eig(first_layer); 
[eig_first, I1] = diag(EL1); 
VL1 = VL1(:,I1); 

[VL2,EL2] = eig(second_layer); 
[eig_second, I2] = diag(EL2);
VL2 = VL2(:,I2); 


C_mat = VL1'*VL2; 
norm_dev(inet1,inet2,idens) = norm(C_mat - eye(N)); 
angle_max(inet1,inet2,idens) = C_mat(end,end); 


% projection of optimal u on the eigenvectors of layer 1
projection11 = zeros(N,N); 
% projection of optimal u on the eigenvectors of layer 2
projection22 = zeros(N,N); 
% projection of optimal u for second layer on the eigenvectors of layer 1
projection12 = zeros(N,N); 


vid_obj1 = VideoWriter(['net1=',num2str(inet1),'net2=',num2str(inet2),'opt_u_dir11_dens=',num2str(idens),'.avi']);
vid_obj2 = VideoWriter(['net1=',num2str(inet1),'net2=',num2str(inet2),'opt_u_dir22_dens=',num2str(idens),'.avi']);
vid_obj3 = VideoWriter(['net1=',num2str(inet1),'net2=',num2str(inet2),'opt_u_dir12_dens=',num2str(idens),'.avi']);


          
 u_layer1 = load(fullfile(dirname,['optimU_L1_inet1=',num2str(inet1),...
                                    '_inet2=',num2str(inet2),'_trial=',num2str(itrial),...
                                    '_dens=',num2str(idens),'.csv']),'delimiter','tab') ;
 u_layer2 = load(fullfile(dirname,['optimU_L2_inet1=',num2str(inet1),...
                                    '_inet2=',num2str(inet2),'_trial=',num2str(itrial),...
                                    '_dens=',num2str(idens),'.csv']),'delimiter','tab') ;
    
 for jj =1:N 
 
     u1_reshape = reshape(u_layer1(jj,:),[N,nt+1]); 
     u2_reshape = reshape(u_layer2(jj,:),[N,nt+1]); 
     
%      u_second = optim_ctrl(:,N+1:2*N); 
     
     % calculate overlap of the control vector with all the eigenvectors 
     % In ideal case, u should be parallel to the eigenmode being excited.
     
     
        projection11(jj,:) = u_first(istep,:)*VL1; 
        projection22(jj,:) = u_second(istep,:)*VL2;
        projection12(jj,:) = u_second(istep,:)*VL1; 
   
        % projection1 and u_first can be transformed into each other by VL1
        % similarly projection2 and u_second. Uncomment the following lines
        % to check this 
        % imagesc(round(projection1*V1' - u_first,15))
        % imagesc(round(projection2*V2' - u_second,15))
end
        
        ener_layer1 = projection11.^2; 
        ener_layer2 = projection22.^2; 
        ener_layer12 = projection12.^2; 
        
        
        ener_layer1 = ener_layer1/max(ener_layer1(:)); 
        ener_layer2 = ener_layer2/max(ener_layer2(:)); 
        ener_layer12 = ener_layer12/max(ener_layer12(:)); 
        
        
        
        
        figure(1)
        fig1 = imagesc(ener_layer1'); 
        colorbar; colormap(summer)
        ylabel('Excited direction relative to layer 1')
        xlabel('Target mode in layer 1')
        set(gca,'FontSize',14)
        axis square
        title(['t = ',num2str(tarray(istep))],'FontWeight','normal')
        drawnow
        open(vid_obj1)
        
        frame1 = getframe(gcf);
        writeVideo(vid_obj1, frame1)
        
        
        figure(2)
       
        fig2 = imagesc(ener_layer2');
        colorbar; colormap(summer)
        set(gca,'FontSize',14)
        ylabel('Excited direction relative to layer 2')
        xlabel('Target mode in layer 2')
        axis square
        title(['t = ',num2str(tarray(istep))],'FontWeight','normal')
        drawnow
        
        frame2 = getframe(gcf);
        open(vid_obj2)
        writeVideo(vid_obj2, frame2)
        
      
        figure(3)
       
        fig3 = imagesc(ener_layer12');
        colorbar; colormap(summer)
        set(gca,'FontSize',14)
        ylabel('Excited direction relative to layer 1')
        xlabel('Target mode in layer 2')
        axis square
        title(['t = ',num2str(tarray(istep))],'FontWeight','normal')
        drawnow
        
        open(vid_obj3) 
        frame3 = getframe(gcf);
        writeVideo(vid_obj3, frame3)
        
        
     
     
     
   
 
end

% close(vid_obj1)
% close(vid_obj2)
% close(vid_obj3)
     
% saveas(gcf, fullfile(fig_dir,['Ener_in_modes_inet1=',num2str(inet1),...
%                                 '_inet2=',num2str(inet2),'_idens=',num2str(idens),...
%                                   '_itrial=',num2str(itrial), '.png'])) ; 











