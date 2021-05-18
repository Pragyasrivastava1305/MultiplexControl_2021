function [state,costate,E] = OneModeSol_funV0(xi, mu,Clm, T, nt, XF)
% function to calculate evolution of weights of state vectors and costate
% vectors along the retained eigen-modes of the two layers. 
% @pragya, Fall 2019

% INPUTS ------------------------------------------------------------------
% xi    :   retained eigenvalue in first layer
% mu    :   retained eigenvalue in second layer
% Clm   :   appropriately scaled dot product of retained eigen-modes
%           (l-th eigenmode of second and m-th eigenmode of first layer )
% T     :   Terminal time
% nt    :   number of intervals on time axis
% XF    :   final state vector in a representation where the two layers are
%           separately diagonalized (2-by-1)in this representation. 

% OUTPUTS -----------------------------------------------------------------
% state :  Dynamical variable of length 2-times-(nt+1). 
%          state(1,:) :  projection of the state of first layer along the
%          retained m-th eigen-mode. 
%          state(2,:) :  projection of the state of first layer along the
%          retained m-th eigen-mode. 
% costate: same as "state", but for costate variables
% E      : energy calculated in one mode approximation


tarray = linspace(0,T,nt+1); 
delt = T/nt; 
state = zeros(2,nt+1); 
costate = zeros(2,nt+1); 


% initial values of projected costate variables ---------------------------
xf1 = XF(1);
xf2 = XF(2); 

if xi ~= mu
    
    den_lam = (exp(2*T*mu)-1) * (exp(2*T*xi)-1) * (mu^2 + xi^2)...
              - 2*(1 + exp(2*T*mu) + exp(2*T*xi) - 4*exp(T*(mu + xi)) + ...
                   exp(2*T*(mu+xi)))*mu*xi; 

    den_sig = den_lam/Clm; 

    lam_factor = 8*exp(T*(mu + xi)) * (mu^2 - xi^2); 

    sig_factor = 4*mu * exp(T*mu) * (mu^2 - xi^2);

    coeff1 = mu * sinh(T*xi) - xi * sinh(T*mu); 

    coeff2 = mu * sinh(T*xi) + xi * (cosh(T*mu) - cosh(T*xi)- sinh(T*mu));

    coeff3 = mu*(exp(2*T*xi) - 1) + xi*(1 + exp(2*T*xi) - 2* exp(T*xi + T*mu));

    coeff4 = (exp(2*T*xi) - 1) * (mu^2 - xi^2); 

    lam0 = -lam_factor * (coeff1 * xf1 + mu * coeff2 * xf2/Clm)/den_lam;

    sig0 = - sig_factor * ( coeff3 * xf1 + coeff4 * xf2/Clm)/den_sig; 

    % Solution for projected state and costate variables ----------------------
    % ------------------------------------------------------------------------

    % projected costate of second layer
    costate(2,:) = exp(-mu*tarray)*sig0; 

    % projected costate of first layer
    costate(1,:) = exp(-xi*tarray)*lam0 + Clm * (exp(-mu*tarray) - exp(-xi*tarray)) ...
                                              * sig0/(mu - xi); 

    % projected state of first layer
    state(1,:) = -sinh(xi*tarray)*lam0/2/xi ...
                 + Clm*( (exp(-mu*tarray) - exp(xi*tarray))/(mu + xi) + ...
                 + sinh(xi*tarray)/xi ) * sig0/2/(mu-xi); 

    % projected state of second layer 
    Factor = ( (exp(-xi*tarray) - exp(mu*tarray))/(mu + xi)  ...
                  +  (exp(xi*tarray) - exp(mu*tarray))/(xi - mu) );


    T1 = - Clm * Factor * lam0/2/xi ;

    T2 = - (Clm^2) * ( (exp(xi*tarray) - exp(mu*tarray))/(xi - mu)  ...
                        - sinh(mu*tarray)/mu ) * sig0/2/(mu^2 - xi^2); 
           
    T3 = (Clm^2) * Factor * sig0/4/xi/(mu - xi);

    state(2,:) = T1 + T2 +T3; 


else
   
    Den_costate = Clm*(1 + 2*(T*xi)^2 - cosh(2*T*xi)); 
    
    costate(2,:) = (8/Den_costate) * ...
                   ( 2*exp(-tarray*xi) *xf2* (xi^3)*sinh(T*xi)/Clm ...
                  - exp(-tarray*xi) *xf1* (xi^2)*( exp(T*xi)*T*xi - sinh(T*xi) ));
              
              
    costate(1,:) = - (8*exp(-tarray*xi) * xi/Den_costate) .* ...
                    (xf2 *xi* (T*xi*cosh(T*xi) + (-1 + 2*tarray*xi - T*xi)) +...
                    xf1 *Clm* (T*xi*(1-tarray*xi)*cosh(T*xi) + (-1 + tarray*xi*(1-T*xi))*sinh(T*xi))...
                    ); 
                    
    Den_state =  (1+ exp(4*T*xi) - 2*exp(2*T*xi) * (1+ 2*(T*xi)^2)); 
    
    state(1,:) = 4*exp(-tarray*xi + T*xi) .* ...
                (xf2* (-1 + exp(2*T*xi))*tarray + T - T*exp(2*tarray*xi)) *xi^2/Den_state/Clm ...
                 ...
                - exp(-tarray*xi + T*xi) *xf1.*...
                (-1 - exp(2*(tarray*xi +T*xi)) + 2*xi*(tarray-T) + exp(2*tarray*xi)*(1+ 2*T*xi) + ...
                exp(2*T*xi)*(1+ 2*tarray*xi*(-1+2*T*xi)))/Den_state;
                

    Term1 = Clm* (tarray*T*xi*cosh(T*xi).*sinh(tarray*xi) - (tarray*T*xi.*cosh(tarray*xi) +(tarray -T).*sinh(tarray*xi)));
    
    Term2 = tarray*xi.*cosh(tarray*xi)*(exp(-T*xi)*T*xi + sinh(T*xi)) + ...
            sinh(tarray*xi) .* (T*xi*(-1 + tarray*xi) *cosh(T*xi) - ...
                                (1+ xi*(tarray-T + tarray*T*xi))*sinh(T*xi)); 
                
                
    state(2,:) = - 4*exp(2*T*xi) *xf1* Term1/Den_state  ...
                 - 4*exp(2*T*xi) *xf2* Term2/Den_state; 
    
end


% calculate energy  and compare with analytic expression ------------------
 E_dens = costate(1,:).^2/4; 
 E = trapz(tarray,E_dens); 
























































end
