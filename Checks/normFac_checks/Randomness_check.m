N=100; 
k = 20; 
G_old = zeros(N,N); 
for itest =1:50
%     rseed(itest) = randi(10e9);
    G2 = net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(k),0.1)); 
    
    Gnew = G2 + G2' - G2.*eye(N);  
    figure; 
    imagesc(Gnew- G_old); colorbar; drawnow
    
    G_old = Gnew;
end