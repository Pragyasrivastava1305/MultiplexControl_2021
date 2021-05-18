function [xbinned, ybinned,min_vec,max_vec] = binned_vectors(xmat, ymat, nbin,lims)

bin_vec = linspace(lims(1),lims(2),nbin+1);

x_master_vec = xmat(:); 
y_master_vec = ymat(:);
% [xbinned, Ibinned] = sort(xmat(:));
% yvec = ymat(:);
% ybinned = yvec(Ibinned); 
% err_vec = zeros(size(yvec,2),1); 

for ibin =1:size(bin_vec,2)-1
    
    % get the indices of xvector falling into different bins
    ind_vec =  (x_master_vec < bin_vec(ibin+1)) & (x_master_vec >= bin_vec(ibin)); 
                        
    % extract values from the original vector in each bin
    if size(find(ind_vec),1) ~=0
        
        % take mean of the points falling inside the interval
        xbinned(ibin) = mean(x_master_vec(ind_vec));
        ybinned(ibin) = median( y_master_vec(ind_vec) ); 
       % err_vec(ibin) = std( y_master_vec(ind_vec) ); 
        min_vec(ibin) = min( y_master_vec(ind_vec) );
        max_vec(ibin) = max( y_master_vec(ind_vec) );
    
    else
        % if no points fall within the interval, take mean of the
        % end-points
        xbinned(ibin) = 0.5*( bin_vec(ibin+1) + bin_vec(ibin) );
        ybinned(ibin) = 0; 
       % err_vec(ibin) = 0; 
        min_vec(ibin) = 0; 
        max_vec(ibin) = 0; 
        
    end
    

end

for ibin =2:size(bin_vec,2)-2
    if ybinned(ibin) ==0
        ybinned(ibin) = 0.5*( ybinned(ibin-1) + ybinned(ibin+1) ); 
        min_vec(ibin) = 0.5*( min_vec(ibin-1) + min_vec(ibin+1) );
        max_vec(ibin) = 0.5*( max_vec(ibin-1) + max_vec(ibin+1) );
        
    end
end

end

