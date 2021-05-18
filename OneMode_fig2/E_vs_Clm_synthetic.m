%close all
m=0;   % to get o
N = 100; 
I = m+1:N-m; 
l = size(I,2);
% inet1 = 2; 
% inet2 = 3; 
for inet1 =1:3
     for inet2 =1:3
      
    cols_n_markers; 
    load(['inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.mat'])
            
    Csubp = reshape(C_mat(I,I),l^2,1);
    Esubp = reshape(E(I,I),l^2,1); 
    plot(abs(Csubp), Esubp,'linestyle','none','Marker',mk,'MarkerFaceColor',col_one,...
                                            'MarkerEdgeColor',col_one,'MarkerSize',12)
    drawnow
    hold on

%     % fastest modes
%     Csub1 = reshape(C_mat(end-m+1:end,end-m+1:end),m^2,1); 
%     Esub1 = reshape(E(end-m+1:end,end-m+1:end),m^2,1); 
%     scatter(abs(Csub1),Esub1,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',0.7)
%     hold on 
%     % slowest modes
%     Csub2 = reshape(C_mat(1:m,1:m),m^2,1); 
%     Esub2 = reshape(E(1:m,1:m),m^2,1); 
%     scatter(abs(Csub2),Esub2,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',0.7)


    set(gca,'XScale','log'); set(gca,'YScale','log')
    grid on ; box on 
    xlabel('|p_i . q_j|'); ylabel('Energy')
    title('Energy vs angle')
     end
end

legend('ER-ER', 'ER-WS', 'ER-BA', 'ER-RG', ...
            'WS-ER', 'WS-WS', 'WS-BA', 'WS-RG', ...
             'BA-ER', 'BA-WS', 'BA-BA', 'BA-RG',...
                'RG-ER', 'RG-WS', 'RG-BA', 'RG-RG' )


