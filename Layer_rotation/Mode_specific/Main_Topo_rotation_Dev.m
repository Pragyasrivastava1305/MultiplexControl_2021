% Code to calculate control channels as a function of layer-alignment
% we store results for each trial
% author:        @pragyasr
% created in:    Fall, 2020
clear all;
% close all;

N = 100;
T = 5;
nt = 50;
tarray = linspace(0,T,nt+1); 
directed = 0; 
stabilize = 0;    % stabilize the duplex if stabilize =1
normalize = 1;    % normalize by max eigen value of first layer when normalize =1
ntrial = 100; 


% resolution of theta axis 
N_rot= 18; 
sv = linspace(0,1,N_rot); 

% clen = 9; 
% blue = [0,0,1];
% green = [0, 1, 0];
% colors_p = [linspace(blue(1),blue(1),clen)', linspace(blue(2),green(2),clen)', linspace(blue(3),green(3),clen)'];



% density of network and the initial state
rho1 = 0.25;
rho2 = 0.25; 
X0 = zeros(2*N,1); 
% assign control input matrix
B = zeros(2*N,2*N); 
B(1:N,1:N) = eye(N); 
itrial =1;



for itrial = 1:ntrial
itrial
% mkdir(['Network_analysis/trial=',num2str(itrial)])
% dirname = ['Network_analysis/trial=',num2str(itrial)]; 
% mkdir(['trial=',num2str(itrial)])
dirname = ['trial=',num2str(itrial)]; 
    
% initiate arrays to store state, optimal input and projection 
x2_array = zeros((nt+1)*2*N, N_rot); 
u2_array = zeros((nt+1)*N, N_rot); 
proj_array =  zeros((nt+1)*N,N_rot); 
mu_array = zeros(N, N_rot); % should not change along N_rot direction


for inet1 =1:4
for inet2 =1:4
% [inet1 inet2]

% build initial duplex
Duplex0 = init_duplex(inet1, inet2, rho1, rho2, N,0);
Duplex = Duplex0; 


% stabilize/normalize as specified
max_eig = max(real(eig(Duplex0))); 
if stabilize ==1
    Duplex = Duplex0/(max_eig +1) -1 ; 
end
% non-dimensionalize time by dividing by the largest eigenvalye of G1
if normalize ==1
    Duplex = Duplex0/max_eig; 
end


% Duplex is the initial unrotated network
% extract layers from the normalized duplex
layer1 = Duplex(1:N,1:N); 
layer2 = Duplex(N+1:2*N, N+1:2*N); 
off_mat = Duplex(N+1:end,1:N); 
% extract scaling of off-diagonal connections  
off_diag = diag(off_mat); 
dscale = off_diag(1);  
  

%------------------- eigen-decomposition of each layer --------------------
[V1,D1] = eig(layer1);      
[D1,I1] = sort(diag(D1));  V1 = V1(:,I1); 
[V2,D2] = eig(layer2); 
[D2,I2] = sort(diag(D2));  V2 = V2(:,I2);
muvec = D2;
xivec = D1; 


%  ROTATION: define the plane of rotation by specifying the axes of plane
ax1  = V1(:,N);     ax2 = V1(:,N-1);
%%%  calculate the components of rotor
rotor = zeros(N,N); 
for i = 1:N
    for j = i+1:N
        rotor(i,j) = ax1(i)*ax2(j) - ax2(i)*ax1(j); 
    end
end
rotor = rotor - rotor';
th = pi/2;

% Define and calculate the rotated Duplex_R
Duplex_R = Duplex; 

alignment = zeros(N,size(sv,2)); 

%%% GENERATE ROTATED NETWORKS AND CALCULATE CONTROL
for js = 1:size(sv,2)
    s = sv(js);
    % define rotation from rotor 
    Rs = expm(-s*th*rotor);
    % rotate the layer
    MatR = Rs*layer2*Rs';
    [VR,DR] = eig(MatR); 
    [DR,IR]= sort(diag(DR));
        
    VR = VR(:,IR);
       
    alignment(:,js) = VR(:,end)'*V1;       
    % define the new duplex 
    Duplex_R(1+N:2*N, 1+N:2*N) = MatR; 
     
    % write matrix for one realization
    if itrial ==1
    writematrix(Duplex_R, fullfile(dirname,['Duplex_inet1=',num2str(inet1),...
                 '_inet2=',num2str(inet2),'_s=',num2str(js), '.csv']),'delimiter','tab') ; 
    end
   
%   imagesc(Duplex_R); colorbar; drawnow
   
    % calculate the state, costate and control input for N-th eigenmode
    XF2 = [zeros(N,1); VR(:,N)]; 
    [x2,u2,v2,n_err2] = min_eng_cont(Duplex, T, B, X0, XF2, nt,0); 
    optim_u2 = u2(:,1:N)'; 
     
    
%   uncomment to check if -VR state has same control energy. It should
%   XF2p = [zeros(N,1); -VR(:,N)]; 
%   [x2p,u2p,v2p,n_err2p] = min_eng_cont(Duplex, T, B, X0, XF2p, nt,0); 
%   optim_u2p = u2p(:,1:N)'; 
         
    for iter = 1:nt+1
        ener_densL2(iter) = optim_u2(:,iter)'*optim_u2(:,iter); 
        
        % calculate the projection of optim_u2 along eigenvectors of the first layer 
        proj(:,iter) = optim_u2(:,iter)'*V1; 
%         ener_densL2p(iter) = optim_u2p(:,iter)'*optim_u2p(:,iter); 
    end

    E2(js) = trapz(tarray,ener_densL2);     
%     E2p(js) = trapz(tarray,ener_densL2p);
    x2_array(:,js) = reshape(x2,[size(x2,1)*size(x2,2),1]); 
    u2_array(:,js) = reshape(optim_u2,[size(optim_u2,1)*size(optim_u2,2),1]);
    proj_array(:,js) = reshape(proj,[size(optim_u2,1)*size(optim_u2,2),1]);
end
 
save(fullfile(dirname,['LayerRotation_inet1=',num2str(inet1),'_inet2=',num2str(inet2),'trial=',num2str(itrial),'_original.mat']))

      
end

end

end


















