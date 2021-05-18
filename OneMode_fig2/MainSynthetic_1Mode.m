% Code to (a) create duplex in different topological classes, (b) change
% the density of second layer to increase spectral distance from first
% layer, and (c) calculate control energy in one-mode limit
% In all these cases, we do parameter sweep on network density 
% and calculate the average controllability as the trace of
% Gramian and difficulty of control from the smallest eigenvalue of
% Gramian. 

% author:       @pragyasr
% created in:    Summer, 2020

clear all;
%close all;
% import python module
%py.importlib.import_module('networkx');
N = 100;
T = 1;
nt =100;
p = 0.1; 
directed = 0; 
ntrial = 10; 

int_vec = 2:2:N-1;  % parameter used to create WS and BA graphs
real_vec = linspace(0.1,1,size(int_vec,2));  % parameter to create ER and RG nets

rho = zeros(size(int_vec,2),4); 
dspec = zeros(size(int_vec,2),4); 
rho_chk = zeros(size(int_vec,2),4); 
lamAmax = zeros(size(int_vec,2),4); 
stdlamA = zeros(size(int_vec,2),4); 
for inet1 =1:3
    for inet2 =1:3
        [inet1 inet2]
% specify the class of networks in first layer and second layer
% 1 == ER, 2 == WS, 3 == BA
% inet1 =1; 
% inet2 =1;

k2 = floor(N/2); 
k4 = floor(N/4); 


if inet1 ==1
    G1 = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N),0.5,pyargs('directed',false))); 
elseif inet1 ==2 
    G1 = net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(k2),0.1)); 
elseif inet1 ==3 
    G1 = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(int16(N),int16(k4)));
end


if inet2 ==1
    loop_vec = real_vec; 
elseif inet2 ==2
    loop_vec = int_vec; 
elseif inet2 ==3
    loop_vec = int_vec; 
end


% make G1 symmetric and traceless
G1 = (G1 + G1') - G1.*eye(N);
% calculate the density of network
rho1 = sum(sum(G1))/N/(N-1); 

B = zeros(2*N,N); 
B(1:N,:) = eye(N); 

Duplex0(1:N,1:N) = G1; 
Duplex0(N+1:2*N,1:N) = eye(N);

clen = 9; 
blue = [0,0,1];
green = [0, 1, 0];
colors_p = [linspace(blue(1),blue(1),clen)', linspace(blue(2),green(2),clen)', linspace(blue(3),green(3),clen)'];
itrial =1;
%for itrial = 1:ntrial
    sd = randperm(1000,1); 
    idens = 24; 
    %for idens = 1:size(loop_vec,2)
    % create second layer with increasing density
        if inet2 ==1
            G2 = net.helper.py_graph2adjmat(py.networkx.erdos_renyi_graph(int16(N),loop_vec(idens),pyargs('directed',false))); 
        elseif inet2 ==2 
            G2 = net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(loop_vec(idens)),0.1,int16(sd))); 
        elseif inet2 ==3 
            G2 = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(int16(N),int16(loop_vec(idens)),int16(sd)));
        end
      
       
       G2 = (G2 + G2') -  G2.*eye(N);  
       Duplex0(N+1:2*N, N+1:2*N) = G2;  
       rho2(idens,itrial) = sum(sum(Duplex0(N+1:2*N,1+N:2*N)))/N/(N-1); 
       
       dspec0(idens,itrial) = spec_dist(G1, G2, 1);
     
       %Duplex = Duplex0; 
       Duplex = Duplex0/eigs(Duplex0(1:N,1:N),1);   
       % extract stabilized layer networks
       layer1 = Duplex(1:N,1:N); 
       layer2 = Duplex(N+1:2*N, N+1:2*N); 
       off_mat = Duplex(N+1:end,1:N); 

       % extract scaling of off-diagonal connections  
       off_diag = diag(off_mat); 
       dscale = off_diag(1);  
       %------------------- eigen-decomposition of each layer --------------------
       [V1,D1] = eig(layer1);      
       [V2,D2] = eig(layer2); 
       [D1,I1] = sort(diag(D1));  V1 = V1(:,I1); 
       [D2,I2] = sort(diag(D2));  V2 = V2(:,I2); 

       xivec = D1; 
       muvec = D2; 
       
       % calculate spectral distance and alignment matrix; 
       dspec1(idens,itrial) = spec_dist(Duplex(1:N,1:N), Duplex(N+1:2*N,N+1:2*N), 1); 
       C_mat = zeros(N,N); 
        for i =1:N
        for j=1:N 
          C_mat(i,j) = V2(:,i)'*V1(:,j); 
        end
        end
    
        dev_mat = V2'*V1 - V1'*V1;
        L2_norm(idens) = norm(dev_mat,2);  
    
        % call function to calculate control energy
        XF = [0,1]; 
   
    for ixi =1:N
        for imu =1:N
            Clm = C_mat(ixi,imu); 
            [state,costate,E(ixi,imu)] =...
                            OneModeSol_funV2(xivec(ixi), muvec(imu),Clm, T, nt, XF); 

        end
    end
    
    save(['inet1=',num2str(inet1),'_inet2=',num2str(inet2),'.mat'])
                %'_rho1=',num2str(round(rho1,1)),'_rho2=',num2str(round(rho2(idens),1)),'.mat'])
    
    end 
end 