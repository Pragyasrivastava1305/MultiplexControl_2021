rho1 = 0.25;
inet1 =2;
N =100;


p1 = get_dens_param(rho1, inet1);  
G1_old = zeros(N,N); 

for isamp =1:50
sd1(isamp) = randperm(5000000,1);
G1_new = net.helper.py_graph2adjmat(py.networkx.barabasi_albert_graph(int16(N),int16(p1)));
G1_new = G1_new +G1_new';

% check density 
rho_chk(isamp) = sum(sum(G1_new))/N/(N-1); 

G1_diff = G1_new - G1_old;

imagesc(G1_new)
title(['trial=',num2str(isamp)])
colorbar
drawnow; 

G1_old = G1_new;

end