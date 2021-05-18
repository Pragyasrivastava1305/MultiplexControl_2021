%clear; close all;  
T= 1; 
Chi = linspace(0,5,50); 
theta = linspace(0,pi/2,91); 
nt = 1000; 
XF = [0 1]';
 
[Xgrid, Ygrid]= meshgrid(cos(theta(1:end-1)), Chi); 


%% with Clm and Chi
for iax = 1:size(Chi,2)
    for iay = 1:size(theta,2)-1
        
        xi = 1; 
        
        mu = Chi(iax)*xi; 
        
        Clm = cos(theta(iay)); 
        %Clm = 1; 
         
        [state,costate,Ean(iax,iay)] = OneModeSol_funV2(xi, mu,Clm, T, nt, XF); 
        
    end
end


surf(Xgrid, Ygrid, Ean)
