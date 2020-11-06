clear all
close all;
ind=25; 
num_trial =1;
ysc = 2;  % if ysc =1, linear scale, if ysc =2, log scale.

% define colors
clen = 4; 
blue = [0,0,1];
green1 = [0, 0.8, 0];
colors_p1 = [linspace(blue(1),blue(1),clen)', linspace(blue(2),green1(2),clen)', linspace(blue(3),green1(3),clen)'];

yellow1 = [204, 153, 0]/255;  yellow2 = [255, 204, 153]/255; 
mauve1  = [204,51,153]/255;    mauve2 = [255,153,255]/255; 
green1 = [0, 153,51]/255;     green2 = [0,255,205]/255;
blue1 = [102,0,255]/255;     blue2 = [0,204,255]/255; 

col1 = [linspace(mauve1(1),mauve2(1),clen)', linspace(mauve1(2),mauve2(2),clen)', linspace(mauve1(3),mauve2(3),clen)'];
col2 = [linspace(green1(1),green2(1),clen)', linspace(green1(2),green2(2),clen)', linspace(green1(3),green2(3),clen)'];
col3 = [linspace(yellow1(1),yellow2(1),clen)', linspace(yellow1(2),yellow2(2),clen)', linspace(yellow1(3),yellow2(3),clen)'];
col4 = [linspace(blue(1),blue2(1),clen)', linspace(blue1(2),blue2(2),clen)', linspace(blue1(3),blue2(3),clen)'];

bcol1 = col4(3,:); 
mcol1 = col1(2,:);
ycol1 = col3(1,:); 
gcol1 = col2(1,:);  

ntop = 4;
strt_ind = 5;
N = 100;
ndens = 25;


for inet1 =1:ntop
%inet2 =1

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
        
    
  
    E1_aver = zeros(N, ndens); 
    E2_aver = zeros(N,ndens);
    EigL1_aver = zeros(N,ndens);
    EigL2_aver = zeros(N,ndens);
    rho2_array = zeros(1,ndens); 
    
    for jtrial =1:num_trial
    dirname = ['trial=',num2str(jtrial),'/inet1=',num2str(inet1),'_inet2=',num2str(inet2)]; 
    
    % load mat file
    load(fullfile(dirname,['topodensity_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'itrial=',num2str(jtrial),'.mat']))
    
    
    E1_aver = E1_aver + E1_array; 
    E2_aver = E2_aver + E2_array; 
    %rho_vec(inet1,inet2) = mean(rho1(inet1,:),2);
    rho_vec(jtrial,inet1) = rho1;
    rho2_array = rho2_array + rho2; 
    
    EigL1_array = zeros(N,ndens); 
    EigL2_array = zeros(N,ndens); 
    
    for idens = 1: size(loop_vec,2)
        SupraA = load(fullfile(dirname,['SupA_inet1=',num2str(inet1),'_inet2=',...
                                             num2str(inet2),'_trial=',num2str(jtrial),'_dens=',num2str(idens),...
                                                                              '.csv']));  

        layer1 = SupraA(1:N,1:N); 
        layer2 = SupraA(N+1:2*N, N+1:2*N); 

        [VL1,DL1] = eig(layer1); 
        [VL2,DL2] = eig(layer2); 
        [EvL1,I1] = sort(diag(DL1));  VL1 = VL1(:,I1); 
        [EvL2,I2] = sort(diag(DL2));  VL2 = VL2(:,I2); 

        EigL1_array(:,idens) = EvL1;
        EigL2_array(:,idens) = EvL2;
    end
    
    EigL1_aver = EigL1_aver + EigL1_array; 
    EigL2_aver = EigL2_aver + EigL2_array; 
    end
     
    % average over trials
    EigL1_aver = EigL1_aver/num_trial; 
    EigL2_aver = EigL2_aver/num_trial; 
    E1_aver = E1_aver/num_trial; 
    E2_aver = E2_aver/num_trial; 
    rho2_array = rho2_array/num_trial; 
    
    % define average control energy and maximum control energy: will have
    % same length as density vec
    averC_L1 = mean(E1_aver); 
    averC_L2 = mean(E2_aver); 
    MaxC_L1 = max(E1_aver); 
    MaxC_L2 = max(E2_aver); 
    
    Evec1 = averC_L1; 
    Evec2 = averC_L2; 
        
        
        figure(1)
        subplot(2,4, inet2)
        plot(rho2_array, Evec1,'o-','color',col_one,'MarkerFaceColor',col_one, 'MarkerEdgeColor',[1 1 1],'linewidth',2,...
                              'MarkerSize',8)
                          
        hold on; 
        set(gca,'fontsize',14,'fontweight','bold')
        %xlabel('\Delta \rho')  
        %ylim([0.83 0.88])
        %ylim([ 70 125])
        xlim([0 1])
        axis square
        drawnow
        
        subplot(2,4, inet2 + ntop)
        plot(rho2_array, Evec2,'o-','color',col_one,'MarkerFaceColor',col_one, 'MarkerEdgeColor',[1 1 1],'linewidth',2,...
                              'MarkerSize',8)
                          
        hold on; 
        set(gca,'fontsize',14,'fontweight','bold')
        %xlabel('\Delta \rho')
        xlim([0 1])
       % ylim([1 1.35])
         %ylim([ 70 122])
        %ylabel('Average control energy for layer 2','FontWeight', 'normal')
        axis square
%       leg  = legend('ER','WS','BA','RG'); 
%       title(leg,'layer 1')
        drawnow
        
        
end 
      
end 


% saveas(gcf,fullfile('figures',['L2_aver&mod_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.fig']))
% saveas(gcf,fullfile('figures',['L2_aver&mod_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.png']))
%          









