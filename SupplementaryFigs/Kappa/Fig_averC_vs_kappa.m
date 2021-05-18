clear all
close all;

N = 100;
num_trial = 5;

kappa_Vec = 1:0.1:10; 
ndens = 11; 
inet1 = 2;
inet2 = 4;
idens = 6; 
        

for jnet1 =1:4
for jnet2 =1:4 

cols_n_markers
E1_aver = zeros(N,size(kappa_Vec,2));
E2_aver = zeros(N,size(kappa_Vec,2));

for jtrial =1:num_trial
  
    load(['topodensity_inet1=',num2str(jnet1),'_inet2=',num2str(jnet2),'itrial=',num2str(jtrial),'.mat'])
    
    E1_aver = E1_aver + E1_array; 
    E2_aver = E2_aver + E2_array; 
    
   
end

E1_aver = E1_aver/num_trial; 
E2_aver = E2_aver/num_trial; 

subplot(1,2,1)
plot(kappa_Vec,mean(E2_aver),'LineStyle','none','marker',mk,'color',col_one,...
    'MarkerFaceColor',col_one, 'MarkerEdgeColor',[1 1 1],'linewidth',1.5,...
                               'MarkerSize',8)
set(gca, 'xscale', 'log')
set(gca, 'yscale','log')
grid on 
drawnow
hold on

subplot(1,2,2)
plot(kappa_Vec,max(E2_aver),'LineStyle','none','marker',mk,'color',col_one,...
    'MarkerFaceColor',col_one, 'MarkerEdgeColor',[1 1 1],'linewidth',1.5,...
                               'MarkerSize',8)

set(gca, 'xscale', 'log')
set(gca, 'yscale','log')
grid on
drawnow
hold on

% for jk =1:size(kappa_Vec,2)
%     if mod(jk,7) ==0
%         plot(muvec, E2_aver(:,jk),'o-');  hold on
%         drawnow
%         
%     end
% end

end
end
subplot(1,2,1)
legend('ER-ER', 'ER-WS', 'ER-BA', 'ER-RG', 'WS-ER', 'WS-WS','WS-BA', 'WS-RG', ...
              'BA-ER', 'BA-WS','BA-BA', 'BA-RG','RG-ER', 'RG-WS','RG-BA', 'RG-RG' )
plot(kappa_Vec, 50*kappa_Vec.^(-2),'-.k' ,'linewidth',2)
set(gca,'fontsize',16)
axis square

subplot(1,2,2)
legend('ER-ER', 'ER-WS', 'ER-BA', 'ER-RG', 'WS-ER', 'WS-WS','WS-BA', 'WS-RG', ...
              'BA-ER', 'BA-WS','BA-BA', 'BA-RG','RG-ER', 'RG-WS','RG-BA', 'RG-RG' )
plot(kappa_Vec, 50*kappa_Vec.^(-2),'-.k' ,'linewidth',2)
set(gca,'fontsize',16)
axis square







