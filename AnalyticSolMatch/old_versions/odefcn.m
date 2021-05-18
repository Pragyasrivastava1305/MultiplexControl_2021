function dydt  = odefcn (t,y,eta,D,C,Dp)
% function to set up dydt to be used as input in odesolver. 

% Inputs: 
% t: independent variable
% y: dependent variable of length 4N. 
%    First N variables : state variables of layer 1
%    Second N variables : state variables of layer 2
%    Third N variables :  costate variables dual to state variables of layer 1
%    Last N variable : costate variables dual to state variables of layer 2

% eta: parameter that fixes significance of state cost relative to input cost

% D : Diagonal matrix containing eigenvalues of layer 1
% Dp: Diagonal matrix containing eigenvalues of layer 2
% C :  Matrix whose (i,j)-th element is dot product of i-th eigenvector of
% layer 1 with j-th eigenvector of layer 2

% Outputs: 
% dydt : a vector of length 4N which contains the RHS of the dynamical
% equations for state and costate variables. 

u = -0.5*y(2*N+1: 3*N); 

dydt(1: N) = D*y(1: N) + u;      % eqn 43 of current version 

dydt(N+1: 2*N) = Dp*y(N+1: 2*N) + C*y(1:N);  % eqn 44 of current version

dydt(2*N+1: 3*N) = - 2*eta*y(1:N) - D*y(2*N+1: 3*N)  - C*y(3*N+1: 4*N); % eqn 47 of current version

dydt(3*N+1: 4*N) = -2*eta*y(N+1: 2*N) - Dp*y(3*N+1: 4*N); 









