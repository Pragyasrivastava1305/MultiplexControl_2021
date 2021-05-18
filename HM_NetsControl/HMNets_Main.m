% check simplex networks
clear all;
close all;

% import colors
color_list; 
col_one = col_list(3,:);

% Set size parameters
N = 100;
ntrial = 40;

N_block = 10;
E = 500;

f_vec= 0.1: 0.05: 0.9; 
g_vec = 2.1: 0.05: 3;

ClustCoef_Av = zeros(size(f_vec,2), size(g_vec,2)); 
DegHet_Av = zeros(size(f_vec,2), size(g_vec,2)); 

for itrial =1:ntrial
for ifrac = 1:size(f_vec,2)
for igam = 1:size(g_vec,2)
     
    A = generate_HMNet(N, N_block, E, f_vec(ifrac), g_vec(igam)); 
 
    deg_vec = sum(A);
    
    rho(ifrac,igam) = sum(sum(A))/N/(N-1); 
    
    ClustCoef_trial(ifrac,igam) = avgClusteringCoefficient(A); 
    
    DegHet_trial(ifrac,igam) = sum(sum(abs(deg_vec'...
                        - deg_vec)))/(N*(N-1))/mean(deg_vec);

end
 
end
 ClustCoef_Av = ClustCoef_trial + ClustCoef_Av; 
 DegHet_Av = DegHet_Av + DegHet_trial; 
 
 figure; 
 imagesc(A); 
 drawnow


end

ClustCoef_Av = ClustCoef_Av/ntrial; 
DegHet_Av = DegHet_Av/ntrial; 

