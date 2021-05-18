close all
jnet1 = 1; 
jnet2 = 2;
jtrial = 10;
averC_vec = [];



for jtrial =1:50
   
    averC_vec = [];

for jnet1 = 1:4; 
    for jnet2 =1:4
dirname = ['trial=',num2str(jtrial)]; 
load(fullfile(dirname,['LayerRotation_inet1=',num2str(jnet1),'_inet2=',num2str(jnet2),'trial=',num2str(jtrial),'_original.mat'])); 

% verify that the alignment of eigenmodes with the rotation 'subspace' does
% not change, and all changes happen only in the rotation plane
% tested for some example inet1, inet2 and imode. In total there will be
% 16*100 figures which can be generated and saved if needed. 

% imode = 50;
for imode =1:N
for jangle = 1:size(sv,2)
    Amat = reshape(alignment(:,jangle),[N,N]);
    
    % alignment of imode-th mode with all the eigenmodes of first layer for
    % all the values of rotation angle.
    mat1(:,jangle) = abs(Amat(imode,:)); 
    
    % angle of imode-th mode with the dominant eigenmode of first layer for
    % all values of jangle
    mat2(imode,jangle) = abs(Amat(imode,end)); 
end

% plot(sv,mat2(imode,:)); hold on; drawnow
% 
end
figure(1)
subplot(5,10,jtrial)
plot(sv, mean(E2_array),'LineWidth',1.5); hold on; 
xlim([0 1]); ylim([60 140])
xlabel('s'); ylabel('average control')
drawnow

% imagesc(align_mat'); 
% figure;
% for jangle =1:18
%     subplot(3,6,jangle)
%     plot(mat2(:,jangle),'LineWidth',2); 
%     xlabel('mode number') 
%     ylabel('alignment with the dominant mode')
%     title(['s = ', num2str(round(sv(jangle),2))])
%     ylim([0 1])
%     drawnow; hold on 
%     
% end

tvec = mean(E2_array); 

averC_vec = [averC_vec tvec(end)]; 
    end
end

figure(2) 
plot(averC_vec, 'o-.', 'LineWidth',1.5); drawnow; hold on


end