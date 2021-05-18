% creates HM nets and calculates their toplogical and spectral properties
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
nbins = 21;

f_vec= 0.1: 0.05: 0.9; 
g_vec = 2.1: 0.05: 3;

clust_coeff_arr = zeros(size(f_vec,2), size(g_vec,2)); 
deg_het_arr = zeros(size(f_vec,2), size(g_vec,2)); 
kurt_arr = zeros(size(f_vec,2), size(g_vec,2)); 
skew_arr = zeros(size(f_vec,2), size(g_vec,2)); 
ex_to_bulk_arr = zeros(size(f_vec,2), size(g_vec,2)); 


eig_array = zeros(N,ntrial); 

for ifrac = 1:size(f_vec,2)
for igam = 1:size(g_vec,2)

    E_hist = zeros(nbins-1,1); 
    xax = zeros(nbins-1,1);
    
    for itrial =1:ntrial
     
    A = generate_HMNet(N, N_block, E, f_vec(ifrac), g_vec(igam)); 
    
    % spectra and properties
    [vv,ee] = eig(A); 
    eig_val = diag(ee)/ee(end,end);
    eig_array(:,itrial) = eig_val; 
    bin_edge = linspace(min(eig_val), max(eig_val),nbins); 
        
    % histogram of eigenvalues
    h = histogram(eig_val,bin_edge,'Visible','off');
    hist_vec = h.Values; 
    hist_vec = hist_vec/sum(hist_vec); 
    E_hist = E_hist + hist_vec'; 
    med_edge = movmean(bin_edge,1); 
    xax = xax + med_edge(2:end)'; 
    
    % proprtiees of eigen-spectrum
    var_vec(itrial) = var(eig_val); 
    kurt_vec(itrial) = kurtosis(eig_val); 
    skew_vec(itrial) = skewness(eig_val); 
   
    % topological properties 
    deg_vec = sum(A);
    rho_vec(itrial) = sum(sum(A))/N/(N-1); 
    clust_coeff(itrial) = avgClusteringCoefficient(A); 
    deg_het(itrial) = sum(sum(abs(deg_vec'...
                        - deg_vec)))/(N*(N-1))/mean(deg_vec);
                    
                    
      
    
    

    end
    
    clust_coeff_arr(ifrac, igam) = mean(clust_coeff);
    deg_het_arr(ifrac,igam) = mean(deg_het); 
    kurt_arr(ifrac,igam) = mean(kurt_vec); 
    skew_arr(ifrac,igam) = mean(skew_vec); 
    var_arr(ifrac,igam) = mean(var_vec); 
    
    E_hist = E_hist/ntrial; 
    
    % define bulk versus extremal values
    ex_to_bulk_arr(ifrac,igam) = sum(E_hist(11:20))/sum(E_hist(1:10)); 
    
    
end
 
 
 figure; 
 imagesc(A); 
 drawnow


end

colormap('bone'); 
subplot(2,3,1) 
imagesc(g_vec, f_vec, deg_het_arr)
axis xy
colorbar
title('Degree Heterogeneity', 'fontweight', 'normal')

subplot(2,3,2) 
imagesc(g_vec, f_vec, clust_coeff_arr)
axis xy
colorbar
title('Cluster Coeff', 'fontweight', 'normal')

subplot(2,3,3) 
imagesc(g_vec,f_vec, ex_to_bulk_arr)
axis xy
colorbar
title('Extreme to bulk', 'fontweight', 'normal')

subplot(2,3,4) 
imagesc(g_vec,f_vec, var_arr)
axis xy
colorbar
title('variance of evs-dist', 'fontweight', 'normal')

subplot(2,3,5) 
imagesc(g_vec,f_vec, skew_arr)
axis xy
colorbar
title('skewness of evs-dist', 'fontweight', 'normal')

subplot(2,3,6) 
imagesc(g_vec,f_vec, kurt_arr)
axis xy
colorbar
title('kurtosis in evs-dist', 'fontweight', 'normal')




















%ClustCoef_Av = ClustCoef_Av/ntrial; 
%DegHet_Av = DegHet_Av/ntrial; 

