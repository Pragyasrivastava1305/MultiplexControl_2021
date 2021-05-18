close all; 
ndens = 25;
def_colors;

% dirname = 'Full_data_mat'; 
fig_opt =1; 
opt_heatmap =0; 
N = 100; 
nt = 50;
% time mean
mean_omega11 = zeros(N,N); 
mean_omega22 = zeros(N,N); 
mean_omega21 = zeros(N,N); 

Teal1 = [178, 223, 219]/255; 
Teal2 = [ 0, 77, 64]/255; 


for jtrial =2:8
    jtrial = 5*jtrial;
%     j_ind
%     jtrial = trial_vec(j_ind);

 t = tiledlayout(4,4);
 t.TileSpacing = 'compact';
for jnet1 = 1:4
for jnet2 = 1:4
%     subplot(4,4, (jnet1-1)*4 + jnet2)
% 
nexttile

for jdens =1:25
    
dirname = ['trial=',num2str(jtrial),'/inet1=',num2str(jnet1),'_inet2=',num2str(jnet2)]; 

SupraA = load(fullfile(dirname,['SupA_inet1=',num2str(jnet1),'_inet2=',...
                                         num2str(jnet2),'_trial=',num2str(jtrial),'_dens=',num2str(jdens),...
                                                                    '.csv']));  
                                                           
% imagesc(SupA) 
% extract layers 
first_layer = SupraA(1:N,1:N); 
second_layer =  SupraA(N+1:2*N, N+1:2*N); 

% eigen-decomposition of the two layers
[VL1,EL1] = eig(first_layer); 
eig_first = sort(diag(EL1)); 


[VL2,EL2] = eig(second_layer); 
eig_second = sort(diag(EL2));

C_mat = VL1'*VL2; 
% norm_dev(inet1,inet2,idens) = norm(C_mat - eye(N)); 
% angle_max(inet1,inet2,idens) = C_mat(end,end); 


% load control inputs for different eigenvectors
optU_L1 = load(fullfile(dirname,['optimU_L1_inet1=',num2str(jnet1),...
                                    '_inet2=',num2str(jnet2),'_trial=',num2str(jtrial),...
                                    '_dens=',num2str(jdens),'.csv'])) ;
                                
optU_L2 = load(fullfile(dirname,['optimU_L2_inet1=',num2str(jnet1),...
                                    '_inet2=',num2str(jnet2),'_trial=',num2str(jtrial),...
                                    '_dens=',num2str(jdens),'.csv'])) ;

% [jnet1 jnet2 jdens]
% projection of optimal u on the eigenvectors of layer 1
omega11 = zeros(nt+1,N); 
% projection of optimal u on the eigenvectors of layer 2
omega22 = zeros(nt+1,N); 
% projection of optimal u for second layer on the eigenvectors of layer 1
omega21 =  zeros(nt+1,N); 

%j_eig = 4;

if opt_heatmap ==1 
    figure((jnet1-1)*4 + jnet2)
    subplot(5,5,jdens)
    imagesc(SupraA); colorbar; drawnow
end
     [jnet1 jnet2 jdens]
for j_eig = 1:N %size(eig_ind,2)
    
     u_first = reshape(optU_L1(j_eig,:),[N,nt+1])';
     u_second = reshape(optU_L2(j_eig,:),[N,nt+1])'; 
     
     % calculate overlap of the control vector with all the eigenvectors 
    
     for istep =1: nt+1
        omega11(istep,:) = u_first(istep,:)*VL1; 
        omega22(istep,:) = u_second(istep,:)*VL2;
        omega21(istep,:) = u_second(istep,:)*VL1; 
   
        % projection1 and u_first can be transformed into each other by VL1
        % similarly projection2 and u_second. Uncomment the following lines
        % to check this 
        % imagesc(round(projection1*V1' - u_first,15))
        % imagesc(round(projection2*V2' - u_second,15))
        
     end
     ener_layer11 = omega11.^2; 
     ener_layer22 = omega22.^2; 
     ener_layer21 = omega21.^2; 
     
     % define mean over time: normalized by the maximum for each 
     mean_omega11(j_eig,:)  =  mean(ener_layer11) / max(mean(ener_layer11));
     mean_omega22(j_eig,:)  =  mean(ener_layer22) / max(mean(ener_layer22));
     mean_omega21(j_eig,:)  =  mean(ener_layer21) / max(mean(ener_layer21));
       
              
  
         
     if fig_opt ==1 && mod(j_eig,20) ==0
          carray = tarr1; 
%         if jnet1 ==1
%             carray = barr1;
%         elseif jnet1 ==2
%             carray = yarr1;
%         elseif jnet1 == 3
%             carray = garr1;
%         else
%             carray = marr1;
%         end
        
        splot = scatter(C_mat(:,j_eig), mean_omega21(j_eig,:),'o','linewidth',2);
        splot.MarkerEdgeColor = [0 0 0]; %[carray(jdens,1), carray(jdens,2),carray(jdens,3)]; 
        splot.MarkerFaceColor = [carray(jdens,1), carray(jdens,2),carray(jdens,3)];
        splot.SizeData = 40;
        splot.LineWidth = 0.5; 
%         splot.Marker = mk_array(jnet2); 
        
        box on; grid on
        hold on; drawnow  
        set(gca,'fontsize',14)
%       xlabel('Alignment');  ylabel('w_{21}^2')
        xlim([-1 1]); ylim([0 1])
      
        % switch off ticklabels for internal plots
        if jnet2 ~= 1;   yticklabels(''); yticks auto; end
        
        if jnet1 ~= 4; xticklabels(''); xticks auto; end
  
        
     end 
                              
        
end

end
if opt_heatmap ==1
saveas(gcf,['SupraA_trial=',num2str(jtrial),'inet1=',num2str(jnet1),'inet2=',num2str(jnet2), '.fig'])
end


% if fig_opt ==1
%     saveas(gcf,['trial=',num2str(jtrial),'inet1=',num2str(jnet1),'inet2=',num2str(jnet2), '.fig'])
% %     close(gcf)
% end

end

end
if fig_opt ==1
    saveas(gcf,['trial=',num2str(jtrial),'.fig'])
    close(gcf)
end

end

