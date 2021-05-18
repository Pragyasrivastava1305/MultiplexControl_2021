function [state,costate,E] = OneModeSol_funV2(xi, mu,Clm, T, nt, XF)
% function to calculate evolution of weights of state vectors and costate
% vectors along the retained eigen-modes of the two layers. 
% @pragya, Fall 2019
% @pragya, Spring 2020: Included more cases
% @pragya: most recent version
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

% Considers four cases differently; (i) xi ~= mu, (ii) xi = mu, 
%                             (iii)xi=0, mu ~= 0, (iv) xi ~=0, mu =0 


tarray = linspace(0,T,nt+1); 
state = zeros(2,nt+1); 
costate = zeros(2,nt+1); 


% initial values of projected costate variables ---------------------------
wf = XF(1);
vf = XF(2);

if (xi ~= mu) && (xi ~= 0) && (mu ~= 0)
    
    fac =  (-1 + exp(2*T*mu))*(-1 + exp(2*T*xi))*(mu^2 + xi^2) ... 
	     - 2*( 1 + exp(2*T*mu) + exp(2*T*xi) - 4*exp(T*(mu + xi)) + exp(2*T*(mu + xi)))*mu*xi; 

    Den_sig = (Clm^2)*fac; 

    Cx_sig = -4*Clm*exp( (T-tarray)*mu )*mu*(mu^2 - xi^2)* ( (-1 + exp(2*T*xi))*mu +...
							     (1 + exp(2*T*xi) - 2*exp( T*(mu+xi)))*xi ...
							   ); 
    Cy_sig =  -4*exp((T-tarray)*mu)*mu*( exp(2*T*xi)-1 )*(mu^2 - xi^2)^2; 

    
    Den_lam = Clm*exp( tarray*(mu + xi) )*fac;    

    Cx_lam = 4*Clm*( mu+xi )*( exp(T*mu + tarray*xi) * mu*(mu-xi) ...
			      + 2*exp( (tarray+T)*mu + 2*T*xi )*mu*xi  ...
		              + 2*exp(2*T*mu + (tarray +T)*xi )*mu*xi  ...
			      - exp(T*xi + tarray*mu) * xi*(mu-xi) ...
			      - exp( tarray*xi + T*(mu + 2*xi) )*mu*(mu+xi)...
			      - exp( tarray*mu + T*(xi + 2*mu) )*xi*(mu+xi)...
			    );

    Cy_lam = 4*mu*( mu^2 - xi^2 )* ( - 2*exp(tarray*mu + T*xi)*xi ...
				     + 2*exp( (tarray+T)*mu + 2*T*xi)*xi ... 
				     + exp( T*mu + tarray*xi )*( mu +xi ) ...		
				     - exp( tarray*xi + T*( mu+2*xi) )*(mu + xi)... 
				   ); 

    Den_X = Den_lam;  

    Cx_X = Clm*( ( 2*exp(T*mu + tarray*xi) - 2*exp( (tarray+T)*mu + 2*tarray*xi) )*mu*(mu-xi)...
                + exp( (tarray + 2*T)*mu + (2*tarray + T)*xi )*((mu - xi)^2) ...
                + 4*exp( 2*T*mu + (tarray+T)*xi )*mu*xi ...
		        + ( 2*exp((tarray+T)*mu + 2*T*xi) - 2*exp(T*mu + tarray*xi + 2*T*xi) )*mu*(mu+xi)...
                + exp( tarray*mu + 2*tarray*xi + T*xi )*( mu^2-xi^2 ) ...
                - exp(tarray*mu + 2*T*mu + T*xi)*( mu+xi )^2 ...
                + exp(tarray*mu + T*xi)*(xi^2-mu^2)...
                );

    Cy_X = -2*mu*(mu^2- xi^2)*( exp( tarray*(mu+ 2*xi) )*(exp(T*mu) - exp(T*xi)) ...
                                + exp( T*mu + tarray*xi )*( -1 + exp(2*T*xi) )...
                                - exp(tarray*mu + T*xi)*( -1 + exp(T*(mu + xi)) )...
                                );


    Den_Y = fac; 

    Cx_Y =  4*Clm*exp(T*(mu + xi)) * ( cosh(T*mu) * ( xi*sinh(tarray*mu) - mu*sinh(tarray*xi) ) + ...
                                       cosh(T*xi) * (-xi*sinh(tarray*mu) + mu*sinh(tarray*xi)) ...
                        + (cosh(tarray*mu) - cosh(tarray*xi))*(-xi*sinh(T*mu) + mu*sinh(T*xi)) ...
                                      ); 

    Cy_Y =  4*exp(T*(mu + xi)) * ( mu*(-cosh(T*mu) + cosh(T*xi) + sinh(T*mu)) * ...
                                   ( -xi* (cosh(tarray*mu) - cosh(tarray*xi) + sinh(tarray*mu)) + ...
                                            mu*sinh(tarray*xi) ) + ... 
                                   ( (mu^2)*cosh(tarray* mu) + (xi^2)*sinh(tarray*mu) - ...
                                      mu*( mu*cosh(tarray*xi) + xi*sinh(tarray*xi) ) ...
                                    )*sinh(T*xi) ...
                                 );
                             
    % Solution for state and costate variables projected along retained eigenmodes ----------------
    % ---------------------------------------------------------------------------------------------
    % projected costate of second layer
    costate(2,:) = (Cx_sig*wf + Cy_sig*vf)./Den_sig;
 
    % projected costate of first layer
    costate(1,:) = (Cx_lam*wf + Cy_lam*vf)./Den_lam; 
 
    % projected state of first layer
    state(1,:) =  (Cx_X*wf + Cy_X*vf)./Den_X; 

    % projected state of second layer 
    state(2,:) = (Cx_Y*wf + Cy_Y*vf)./Den_Y;
   
                                  
elseif (xi == mu) && (xi ~= 0) && (mu ~= 0)
   
   fac =  (1 + 2*(T*xi)^2 - cosh(2*T*xi));

   
   Den_sig = (Clm^2)*exp(tarray*xi)*fac; 

   Cx_sig = -8*Clm*(xi^2)*( exp(T*xi)*T*xi - sinh(T*xi) );
   
   Cy_sig = 16*(xi^3)*sinh(T*xi); 

   
   Den_lam = Clm*exp( (tarray +T)*xi )*fac;    

   Cx_lam = 4*Clm*xi*( -1 + tarray*xi -T*xi + exp(2*T*xi)*( 1-(tarray+T)*xi + 2*tarray*T*(xi^2) ) ); 

   Cy_lam = -4*(xi^2)*( 1 - 2*tarray*xi + 2*T*xi + exp(2*T*xi)*( -1 + 2*tarray*xi ) ); 

    
   Den_X = Clm*exp( T*xi )*fac;  

   Cx_X = Clm*( tarray*xi*( 1 + exp(2*T*xi)*( -1 + 2*T*xi) ).*cosh(tarray*xi) ... 
               + (1 - tarray*xi + 2*T*xi + exp(2*T*xi).*( -1 + tarray*xi*(1 - 2*T*xi) ) ).*sinh(tarray*xi) ...
               ); 

   Cy_X = 2*exp( -tarray*xi).* ( tarray - exp(2*T*xi)*tarray + ( -1 + exp(2*tarray*xi) )*T )*(xi^2);

   
   Den_Y = 2*exp( (tarray+T)*xi )*fac; 

   Cx_Y = Clm*( ( -1 + exp(2*tarray*xi) )*( -1 + exp(2*T*xi ) )*T ...
           + tarray.*( -1 - exp( 2*(tarray + T)*xi ) + exp( 2*T*xi)*(1 - 2*T*xi) + exp( 2*tarray*xi)*(1 + 2*T*xi) )...
              );

   Cy_Y =  -1 - exp( 2*(tarray + T)*xi ) - 2*tarray*xi + 2*T*xi + exp( 2*T*xi )*( 1 + 2*tarray*xi ) ...
             + exp( 2*tarray*xi ).*( 1 + 2*T*xi*( -1 + 2*tarray*xi ) );
                             
   % Solution for state and costate variables projected along retained eigenmodes ----------------
   % ---------------------------------------------------------------------------------------------
   % projected costate of second layer
   costate(2,:) = (Cx_sig*wf + Cy_sig*vf)./Den_sig;
 
   % projected costate of first layer
   costate(1,:) = (Cx_lam*wf + Cy_lam*vf)./Den_lam; 
 
   % projected state of first layer
   state(1,:) =  (Cx_X*wf + Cy_X*vf)./Den_X; 

   % projected state of second layer 
   state(2,:) = (Cx_Y*wf + Cy_Y*vf)./Den_Y;
   
   
   
elseif (xi ~= mu) && (xi ==0) && (mu ~=0) 
    
   fac =  ( -2 + T*mu*coth(0.5*T*mu) );

   
   Den_sig = (Clm^2)*exp(tarray*mu)*fac; 

   Cx_sig = -Clm*(mu^2)*( 1 - exp(T*mu)+T*mu )*( csch(0.5*T) )^2; 
   
   Cy_sig = -T*(mu^4)*( csch( 0.5*T*mu) )^2; 

   
   Den_lam = Clm*exp( tarray*mu )*fac;    

   Cx_lam = -Clm*mu*( -2 + 2*exp( tarray*mu ) - 2*coth(0.5*T*mu) + T*mu*(csch(0.5*T*mu))^2 );
   
   Cy_lam = (mu^2)*( 2*exp(tarray*mu)*( -1 + coth(0.5*T*mu) ) - T*mu*( csch(0.5*T*mu) )^2 ); 

    
   Den_X = 2*Clm*exp( tarray*mu )*fac;  

   Cx_X = Clm*( 2 + 2*exp(tarray*mu).*( -1 + tarray*mu ) - 2*( -1 + exp(tarray*mu) )*coth(0.5*T*mu)...
                + (-1 + exp(tarray*mu) )*T*mu*( csch(0.5*T*mu) )^2 );

   Cy_X = -(mu^2)*( -2*exp(tarray*mu).*tarray + 2*exp( tarray*mu ).*tarray*coth(0.5*T*mu)...
                    - (-1 + exp(tarray*mu) )*T*( csch(0.5*T^2) )^2 ); 

   
   Den_Y = 2*mu*fac; 

   Cx_Y = Clm*( (csch(0.5*T*mu))^2 )*( tarray*mu - T*mu + T*mu*cosh(tarray*mu) ...
                 - tarray*mu*cosh( T*mu) - sinh(tarray*mu) + sinh((tarray - T)*mu) + sinh(T*mu)...
                                      );
   
   Cy_Y =  mu*( (csch(0.5*T*mu))^2 )*( 1 + tarray*mu - T*mu + ( -1 + T*mu )*cosh(tarray*mu) ...
                  + cosh( (tarray - T)*mu ) - (1 + tarray*mu)*cosh(T*mu) - sinh( tarray*mu )...
                  + sinh( (tarray - T)*mu ) + ( 1 + tarray*mu )*sinh(T*mu) ...
                                      ); 
   
   % Solution for state and costate variables projected along retained eigenmodes ----------------
   % ---------------------------------------------------------------------------------------------
   % projected costate of second layer
   costate(2,:) = (Cx_sig*wf + Cy_sig*vf)./Den_sig;
 
   % projected costate of first layer
   costate(1,:) = (Cx_lam*wf + Cy_lam*vf)./Den_lam; 
 
   % projected state of first layer
   state(1,:) =  (Cx_X*wf + Cy_X*vf)./Den_X; 

   % projected state of second layer 
   state(2,:) = (Cx_Y*wf + Cy_Y*vf)./Den_Y;
    
     
 elseif (xi ~= mu) && (xi ~=0) && (mu ==0) 
       
   fac =  (2 - 2*cosh(T*xi) + T*xi*sinh(T*xi)); 
   
   Den_sig = (Clm^2)*fac;  

   Cx_sig = 2*Clm*(xi^2)*( -1 + cosh(T*xi) ); 
   
   Cy_sig = -2*(xi^3)*sinh(T*xi); 

   
   Den_lam = Clm*exp( tarray*xi )*fac;    

   Cx_lam = -2*Clm*xi*( 1 - exp(tarray*xi) - exp(T*xi) +T*xi + exp(tarray*xi)*cosh(T*xi) );
   
   Cy_lam = 2*xi*( xi - exp(T*xi)*xi + exp(tarray*xi)*xi*sinh(T*xi) ); 

    
   Den_X = -Clm*fac;  

   Cx_X = Clm*( -1 + cosh( tarray*xi ) - cosh( (tarray - T)*xi ) + cosh(T*xi) - T*xi*sinh(tarray*xi) ); 
   
   Cy_X = xi*( sinh(tarray*xi) - sinh( (tarray - T)*xi ) - sinh(T*xi) ); 

   
   Den_Y = xi*fac; 

   Cx_Y = Clm*( tarray*xi - T*xi + T*xi*cosh( tarray*xi ) - tarray*xi*cosh(T*xi) - sinh(tarray*xi) ...
                + sinh((tarray - T)*xi) + sinh(T*xi) ); 
   
   Cy_Y =  xi*( 1 - cosh(tarray*xi) + cosh((tarray - T)*xi) - cosh(T*xi) + tarray*xi*sinh(T*xi)); 
   
   % Solution for state and costate variables projected along retained eigenmodes ----------------
   % ---------------------------------------------------------------------------------------------
   % projected costate of second layer
   costate(2,:) = (Cx_sig*wf + Cy_sig*vf)./Den_sig;
 
   % projected costate of first layer
   costate(1,:) = (Cx_lam*wf + Cy_lam*vf)./Den_lam; 
 
   % projected state of first layer
   state(1,:) =  (Cx_X*wf + Cy_X*vf)./Den_X; 

   % projected state of second layer 
   state(2,:) = (Cx_Y*wf + Cy_Y*vf)./Den_Y;
    
   
    
end


% calculate energy  and compare with analytic expression ------------------
 E_dens = costate(1,:).^2/4; 
 E = trapz(tarray,E_dens); 
























































end
