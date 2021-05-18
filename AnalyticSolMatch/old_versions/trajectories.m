%function [x1, x2, lam,sig] = trajectories(A1,A2,kappa,L0,XF,T,nt) 
% function to calculate trajectories according to analytical 
% expressions (update reference to the notes)
% created by @pragya in Fall 2019 
% L0 : [lam0; sig0];

% Inputs ======================================================
% A1	: adjacency matrix of first layer, stabilized by the largest eigenvalue of 
%	  duplex (N-by-N)
% A2    : same as A1 but for second layer (N-ny-N)
% Both A1 and A2 are symmetric and stabilized matrices, thus all eigenvalues 
% are real and negative.
% kappa : interlayer connectivity matrix between two layers (N-by-N)
% X0	: initial state vector (2N-by-1), first N -> nodes in layer 1
%	  second N -> nodes in layer 2
% XF 	: Terminal state vector
% T 	: time horizon 
% nt    : number of time points -1
% rd	: rounding decimal place

% Outputs =====================================================
% x1	: N-by-(nt+1) array containing state trajectories of nodes in layer 1
% x2 	: N-by-(nt+1) array containing state trajectories of nodes in layer 2	
% lam   : N-by-(nt+1) array of costate trajectories dual to nodes in layer 1
% mu    : N-by-(nt+1) array of costate trajectories dual to nodes in layer 2

A1 = Sym_WS;
A2 = Sym_BA;
N = size(A1,1);
tarray = linspace(0,T,nt+1); 
sig = zeros(N,nt+1);
lam = zeros(N,nt+1);
x1 = zeros(N,nt+1);
x2 = zeros(N,nt+1);

[V1,D1] = eig(A1); 
[V2,D2] = eig(A2); 

% xivec and muvec are naturally ordered in ascending order, thus most negative
% eigenvalue (largest in magnitude) will be the first one
xivec = spdiags(D1); 
muvec = spdiags(D2); 

% obtain weights of initial and final state vectors
for ivec =1:N
% 	x10(ivec) = dot(X0(1:N),V1(:,ivec)); 
% 	x20(ivec) = dot(X0(N+1:2*N),V2(:,ivec)); 
    lam0(ivec) = dot(L0(1:N),V1(:,ivec)); 
 	sig0(ivec) = dot(L0(N+1:2*N),V2(:,ivec)); 
	x1f(ivec) = dot(XF(1:N),V1(:,ivec)); 
	x2f(ivec) = dot(XF(N+1:2*N),V2(:,ivec));
end
lam0 = lam0';
sig0 = sig0';
% define arrays to store the trajectors of weights
sig_inter = zeros(N,nt+1);
lam_inter = zeros(N,nt+1);
x1_inter = zeros(N,nt+1);
x2_inter = zeros(N,nt+1);

C = zeros(N,N); 

% calculate the matrix C 
for irow=1:N
	for jcol=1:N
	C(irow,jcol) = V2(:,irow)'*kappa*V1(:,jcol); 
	end
end

% calculate matrices to compute state trajectory of the second layer nodes
M = zeros(N,N,nt+1); 
K1 = zeros(N,N,nt+1); 
K2 = zeros(N,N,nt+1);

% calculate matrices appearing in the solution x2
for ii=1:N
	for jj=1:N
		for iT =1:nt+1
        		M(ii,jj,iT) = -C(ii,jj)*exp(muvec(ii)*tarray(iT))*( (exp(-(xivec(jj) + muvec(ii))*tarray(iT)) -1) / (xivec(jj) + muvec(ii))...
				  + (exp((xivec(jj) - muvec(ii))*tarray(iT)) -1) / (xivec(jj) - muvec(ii)) ); 
        		
			K1(ii,jj,iT) = C(ii,jj)*C(jj,jj)*exp(muvec(ii)*tarray(iT))*( (exp(-(muvec(jj) + muvec(ii))*tarray(iT)) -1) / (-muvec(jj) - muvec(ii))...
				  + (exp((xivec(jj) - muvec(ii))*tarray(iT)) -1) / (xivec(jj) - muvec(ii)) )/2/(muvec(jj)^2-xivec(jj)^2); 
			
			K2(ii,jj,iT) = -C(ii,jj)*C(jj,jj)*exp(muvec(ii)*tarray(iT))*( exp((xivec(jj) - muvec(ii))*tarray(iT)) ...
					* (tarray(iT)*(xivec(jj)-muvec(ii)) -1) +1 )/2/(muvec(jj) - xivec(jj)) / ((xivec(jj)- muvec(ii))^2); 
                
		end
	end
end

%MT = M(:,:,end);
%KT = K1(:,:,end) + K2(:,:,end); 

% Calcualte the initial values of costate variables
%[sig0, lam0]= initial_costate(x1f,x2f,xivec,muvec,C,MT,KT,N,T); 

% Trajectories of weights
for im =1:N
	sig_inter(im,:) = sig0(im)*exp(-muvec(im)*tarray); 

	lam_inter(im,:) = lam0(im)*exp(-xivec(im)*tarray) ...
			    + C(im,im)*sig0(im)*(exp(-muvec(im)*tarray) ...
					    + exp(xivec(im)*tarray))/(muvec(im) - xivec(im)) ;

	x1_inter(im,:) = lam0(im)*(exp(-xivec(im)*tarray) - exp(xivec(im)*tarray))/4/xivec(im)...
			   +C(im, im)*sig0(im)*(exp(-muvec(im)*tarray) -exp(xivec(im)*tarray))/2/(muvec(im)^2 - xivec(im)^2);
			    
end

% state trajectory of second layer nodes
lvec  = lam0./xivec/4; 

for iT =1:nt+1
	x2_inter(:,iT) = M(:,:,iT)*lvec + (K1(:,:,iT) + K2(:,:,iT))*sig0; 
end


% convert all the weight functions to state varibles 
for iT =1:nt+1
    x1(:,iT) = V1*x1_inter(:,iT); 
    x2(:,iT) = V2*x2_inter(:,iT);
    lam(:,iT) = V1*lam_inter(:,iT); 
    sig(:,iT) = V2*sig_inter(:,iT); 
end



% function [sig0,lam0] = initial_costate(x1f,x2f,xivec,muvec,C,MT,KT,N,T)
% % calculate initial values of costate variables from here
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g1 = zeros(N,1); 
% g2 = zeros(N,1);
% 
% % solution for sig0
% for ic =1:N 
% 	fac = 4*xivec(ic)/(exp(-xivec(ic)*T) - exp(xivec(ic)*T)); 
% 	g1(ic) = fac*x1f(ic); 
% 
% 	g2(ic) = fac*(C(ic,ic)*T*exp(xivec(ic)*T)/(muvec(ic) - xivec(ic))/2  ...
% 		          - C(ic,ic)*( exp(-muvec(ic)*T) - exp(xivec(ic)*T) )/(muvec(ic)^2 - xivec(ic)^2)/2 );     
% end
% 
% Vtilde = MT*(g1./xivec/4); 
% mat =zeros(N,N); 
% 
% for i=1:N
% 	for j=1:N
% 		mat(i,j) = MT(i,j)*g2(j)/4/xivec(j); 
% 	end 	
% end
% 
% Mtilde = mat + KT;
% sig0 = pinv(Mtilde)*(x2f' - Vtilde); 
% lam0 = g1 + g2.*sig0; 
% 
% 






































 
