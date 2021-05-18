function [param_gen] = get_dens_param(req_dens, net_num)
% This function outputs the parameter to generate the network of 
% required density specified in req_dens and topology specified in net_num 


data_matrix = load('rho_vs_netparams.csv'); 

x_axis = data_matrix(:,net_num); 
y_axis = data_matrix(:,net_num + 4); 

if (net_num ~=3)
    param_gen = interp1(y_axis,x_axis,req_dens); 

else 
    fun = @(x) -0.000202*x.^2 + 0.0202*x- req_dens;
    param_gen = fzero(fun, 10); 
    
    % fun^{-1} is a double-valued function of the dependent variable. 
    % param_gen gives the lower root. 
 
end

if net_num ==2 || net_num ==3
    param_gen = floor(param_gen); 
end



end

