% Code to identify input eigenmodes that carry optimal energy of 
% control, when the target state is a combination of eigenmodes of the
% second layer

clear all
close all

py.importlib.import_module('networkx');

% Network size 
N = 100; 
m = 10;

% Time parameters
T = 5; 
nt = 50; 
tarray =  linspace(0,T,nt+1); 
Ntrial = 100; 
B = zeros(2*N, N); 
B(1:N,1:N) =  eye(N); 

mov = VideoWriter('U_ApproxSys_byMode10.avi'); 
% create original duplex by calling the following script
rho1 = 0.25; 
rho2 = rho1;
% col = [184, 115, 51]/255;
color = [0.8, 0.8, 0.8]; 

for inet1 =1:4
 
Duplex0 = init_duplex_identical_layers(inet1,rho1,N, 0); 

scale_fac = eigs(Duplex0(1:N, 1:N),1 ); 

Duplex = Duplex0/scale_fac; 

L1 = Duplex(1:N, 1:N); 
% set the two layers to be the same
L2 = Duplex(1+N:2*N, 1+N:2*N); 

[V1,D1] = eig(L1);      
% [V2,D2] = eig(L2); 
 
[D1,I1] = sort(diag(D1));  V1 = V1(:,I1); 
%[D2,I2] = sort(diag(D2));  V2 = V2(:,I2); 

P = V1; 
Q = P;

% Set initial conditions
X1_init = zeros(N,1);    
X2_init = zeros(N,1);
X0 = [X1_init; X2_init];

u_overlap = zeros(nt+1,N);  % collect overlap of u with each eigenvector in Q. 

% pick the eigenvectors that contribute to final state

m2(:,inet1) = sort(randi(N,[1,m]));

vec_P = zeros(N,1);     
vec_Q = zeros(N,1);
rvec = 0.5*rand(m,1) + 0.5;

for col = 1:m
    vec_Q = vec_Q +  rvec(col)*Q(:,m2(col,inet1)); 
end

vec_Q = vec_Q/norm(vec_Q); 

XF = [vec_P; vec_Q]; 
XF = XF/norm(XF);

[state,u,v,n_err] = min_eng_cont(Duplex, T, B, X0, XF, nt,0); 

for iter = 1:nt+1
    
    u_overlap(iter,:) = u(iter,:)*P; 
%      subplot(2,2,inet1)
    u_sq(iter,:) = u_overlap(iter,:).^2;
%     plot(1:N, u_sq(iter,:),'.'); 
%     hold on;
%     drawnow
    
end
    u_sq_norm = u_sq/max(u_sq(:)); 
    subplot(2,2,inet1)
       
%   imagesc(tarray,1:N,u_sq_norm');
   
    for iter =1:nt+1
    plot(1:N, u_sq_norm(iter,:),'sq','MarkerEdgeColor','k','MarkerFaceColor',color,...
        'LineWidth',1); hold on
    xlabel('time','FontSize',18)
    ylabel('mode number','FontSize',18);
    set(gca,'FontSize',16)
%     colorbar
%     caxis([  0 1])
%     colormap('hot'); 
    drawnow
    end
    
    plot(m2(:, inet1), 0*ones(m,1),'o', 'MarkerSize',10, 'MarkerfaceColor'...
                ,'g', 'MarkerEdgeColor','k');  hold on


end













