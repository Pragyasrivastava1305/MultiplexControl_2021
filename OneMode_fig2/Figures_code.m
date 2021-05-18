%close all
figure(2)
subplot(2,3,1)
% plot original
imagesc(Duplex0);  colorbar; 
set(gca,'FontSize',14); set(gca,'FontSize',12)
xlabel('Brain regions'); ylabel('Brain regions')
title('Original duplex')
m= 5;
% plot stabilized duplex
subplot(2,3,2)
imagesc(Duplex); colorbar; 
set(gca,'FontSize',14); set(gca,'FontSize',12)
xlabel('Brain regions'); ylabel('Brain regions')
title('Stabilized duplex')


% plot eigen-value distributions 
subplot(2,3,3) 
histogram(xivec,'FaceAlpha',0.25);
hold on
histogram(muvec,'FaceAlpha',0.25);
legend('layer1','layer2')
title('histogram of eigen-values')


% plot energy map
subplot(2,3,4) 
imagesc(E); colorbar; 
set(gca,'FontSize',14); set(gca,'FontSize',12)
xlabel('Brain regions'); ylabel('Brain regions')
title('modal map of energy')

% plot angle map
subplot(2,3,5) 
imagesc(C_mat); colorbar
set(gca,'FontSize',14); set(gca,'FontSize',12)
xlabel('Brain regions'); ylabel('Brain regions')
title('Projection of eigen-vectors')

% plot energy against angle
subplot(2,3,6) 
%%
% C_vec = reshape(C_mat,[N^2,1]); 
% E_vec = reshape(E,[N^2,1]); 
% scatter(abs(C_vec),E_vec,'MarkerFaceColor',c2,'MarkerEdgeColor',c2,'MarkerFaceAlpha',0.3)
% hold on;
% bulk modes
m=0;
I = m+1:N-m; 
l = size(I,2);
Csubp = reshape(C_mat(I,I),l^2,1);
Esubp = reshape(E(I,I),l^2,1); 
scatter(abs(Csubp), Esubp,'.')
hold on

% fastest modes
Csub1 = reshape(C_mat(end-m+1:end,end-m+1:end),m^2,1); 
Esub1 = reshape(E(end-m+1:end,end-m+1:end),m^2,1); 
scatter(abs(Csub1),Esub1,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',0.7)
hold on 
% slowest modes
Csub2 = reshape(C_mat(1:m,1:m),m^2,1); 
Esub2 = reshape(E(1:m,1:m),m^2,1); 
scatter(abs(Csub2),Esub2,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',0.7)


set(gca,'XScale','linear'); set(gca,'YScale','log')
grid on ; box on 
xlabel('|p_i . q_j|'); ylabel('Energy')
title('Energy vs angle')
%% 
% figure(2)
% subplot(1,2,1)
% violin(log(E_vec))
% axis square
% xticks([])
% ylabel ('Energy values (log scale)')
% subplot(1,2,2)
% evs = eig(W)
% violin(log(evs))
% axis square 
% ylabel ('Eigenvalues of Gramian (log scale)')
