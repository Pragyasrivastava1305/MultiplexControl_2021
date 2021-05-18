clear all
jtrial = 5;
jnet1 = 1;
jnet2 = 1;
jdens = 25;
N =100; 
nt = 50;
dirname = ['trial=',num2str(jtrial),'/inet1=',num2str(jnet1),'_inet2=',num2str(jnet2)]; 

% load the supradjacencey matrix
SupraA = load(fullfile(dirname,['SupA_inet1=',num2str(jnet1),'_inet2=',...
                                         num2str(jnet2),'_trial=',num2str(jtrial),'_dens=',num2str(jdens),...
                                                                    '.csv']));  
% extract layers and calculate eigenvalues/vectors
first_layer = SupraA(1:N,1:N); 
second_layer =  SupraA(N+1:2*N, N+1:2*N); 

% eigen-decomposition of the two layers
[VL1,EL1] = eig(first_layer); 
[eig_first, I1] = sort(diag(EL1)); 
VL1 = VL1(:,I1); 

[VL2,EL2] = eig(second_layer); 
[eig_second, I2] = sort(diag(EL2));
VL2 = VL2(:,I2); 

C_mat = VL1'*VL2; 
norm_dev(jnet1,jnet2,jdens) = norm(C_mat - eye(N)); 
angle_max(jnet1,jnet2,jdens) = C_mat(end,end); 
                                                       

% load data of optimal control signal
optU_L1 = load(fullfile(dirname,['optimU_L1_inet1=',num2str(jnet1),...
                                    '_inet2=',num2str(jnet2),'_trial=',num2str(jtrial),...
                                    '_dens=',num2str(jdens),'.csv'])) ;
                                
optU_L2 = load(fullfile(dirname,['optimU_L2_inet1=',num2str(jnet1),...
                                    '_inet2=',num2str(jnet2),'_trial=',num2str(jtrial),...
                                    '_dens=',num2str(jdens),'.csv'])) ;
                                
% time mean
mean_omega11 = zeros(N,N); 
mean_omega22 = zeros(N,N); 
mean_omega21 = zeros(N,N); 

% initialize movie 
vid_obj1 = VideoWriter('test.avi');

for itarget = 1:N
% specify target mode-index 
% itarget = 100; 
 
u1_reshape = reshape(optU_L1(itarget,:),[N,nt+1]); 
u2_reshape = reshape(optU_L2(itarget,:),[N,nt+1]); 
     
% calculate omega as defined in the paper 
     
omega11 = VL1'*u1_reshape; 
omega22 = VL2'*u2_reshape;
omega21 = VL1'*u2_reshape; 
   
        % check if transformation is done right
%         subplot(1,2,1)
%         imagesc(round(VL1*omega11 - u1_reshape,12))
%         colorbar
%         subplot(1,2,2)
%         imagesc(round(VL2*omega22 - u2_reshape,12))
%         colorbar
%         drawnow
        ener_layer11 = omega11.^2; 
        ener_layer22 = omega22.^2; 
        ener_layer21 = omega21.^2; 
            
       mean_omega11(itarget,:)  =  mean(ener_layer11,2)'/ max(mean(ener_layer11,2));
       mean_omega22(itarget,:)  =  mean(ener_layer22,2)'/ max(mean(ener_layer22,2));
       mean_omega21(itarget,:)  =  mean(ener_layer21,2)'/ max(mean(ener_layer21,2));
       
              
        
        figure(1)
        fig1 = imagesc(ener_layer11'); 
        colorbar; colormap(summer)
%         ylabel('Excited direction relative to layer 1')
%         xlabel('Target mode in layer 1')
%         set(gca,'FontSize',14)
%         axis square
%         title(['t = ',num2str(tarray(istep))],'FontWeight','normal')
%         drawnow
%         open(vid_obj1)
%         
%         frame1 = getframe(gcf);
%         writeVideo(vid_obj1, frame1)
        
        
        figure(2)
       
        fig2 = imagesc(ener_layer22');
        colorbar; colormap(summer)
        set(gca,'FontSize',14)
%         ylabel('Excited direction relative to layer 2')
%         xlabel('Target mode in layer 2')
%         axis square
%         title(['t = ',num2str(tarray(istep))],'FontWeight','normal')
%         drawnow
%         
%         frame2 = getframe(gcf);
%         open(vid_obj2)
%         writeVideo(vid_obj2, frame2)
        
      
        figure(3)
       
        fig3 = imagesc(ener_layer21');
        colorbar; colormap(summer)
        set(gca,'FontSize',14)
%         ylabel('Excited direction relative to layer 1')
%         xlabel('Target mode in layer 2')
%         axis square
%         title(['t = ',num2str(tarray(istep))],'FontWeight','normal')
%         drawnow
        
%         open(vid_obj3) 
%         frame3 = getframe(gcf);
%         writeVideo(vid_obj3, frame3)
%         
                                                                            
                                                                            
  end                                                                         
                                                                            
                                                                            
                                                                            
                                                                            
                                    