% This code is to check if the method of direct analytical solution 
% works for simple case: simplex with all node input. 
clear all;

% import python module
py.importlib.import_module('networkx');

N = 100;

p = 0.1; 

directed = 0; 

% Create duplex : layer 1 is WS and layer 2 is BA of comparable density (0.5)
G_WS =net.helper.py_graph2adjmat(py.networkx.watts_strogatz_graph(int16(N),int16(51),0.1));
Sym_WS = 0.5*(G_WS + G_WS');

% Stabilize the network 
Sym_WS = Sym_WS/(eigs(Sym_WS,1)+1) - eye(size(Sym_WS));

% eigenvalues and vectors of first and second layer networks
[V1,D1] = eig(Sym_WS); 
[D1,I1] = sort(diag(D1));  
V1 = V1(:,I1); 

nT = 1000; 
% arrays to store solutions 
x1sol = zeros(N,nT+1); 
lamsol= zeros(N,nT+1); 

T = 1; 
tarray = linspace(0,T,nT+1);
delt = T/nT; 
B = eye(N); 

for iter =1:100
    % initialize arrays for analytical solutions
    x1sol = zeros(N,nT+1); 
    lamsol = zeros(N,nT+1);
        
    % initial state at origin
    X0 = zeros(N,1);
    % generate random final states
    rng(iter)
    XF = randn(N,1);
    XF = XF/norm(XF);
    
    % get numerical solution
    [x,u,v,n_err] = min_eng_cont(Sym_WS, T, B, X0, XF,nT, 0); 
    lam_num = v(N+1:end,:); 
    
    % get initial and final weights along eigenvectors
    x0 = pinv(V1)*X0; 
    xf = pinv(V1)*XF;
    
    % initial state of costate variables
    lam0 = -2*xf.*D1./(sinh(D1*T));   
    
    for iT =1:nT+1
        x1sol(:,iT) = xf.*sinh(D1*tarray(iT))./sinh(D1*T); 
        lamsol(:,iT) = lam0.*exp(-D1.*tarray(iT));     
        x1comp(:,iT) = V1*x1sol(:,iT); 
        lamcomp(:,iT) = V1*lamsol(:,iT); 
    end

    % check that exp(-Sym_BA*T)*lc0 gives correct result as well 
    lc0 = V1*lam0;
    lam_chk = zeros(N,nT+1); 
    lam_chk(:,1) = lc0; 
    for it =1:nT
    lam_chk(:,it+1) = expm(-Sym_WS*T/nT)*lam_chk(:,it);
    end

    %imagesc(lam_chk- lamcomp);colorbar
    [iter max(max(abs(lam_num-lam_chk))) max(max(abs(x'-x1comp))) ]
    
    for imode = 1:N
        Enum(iter,imode) = trapz(tarray,u(:,imode).^2);         
        Ean(iter,imode) = (lam0(imode)^2)*(1- exp(-2*D1(imode)*T))/8/D1(imode); 
    end

    EN(iter) = sum(Enum(iter,:)); 
    EA(iter) = sum(Ean(iter,:)); 
end

%%%% energy calculated like this matches with the numerics. 

