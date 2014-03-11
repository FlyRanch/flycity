function [ w_avg, w_var, w_std ] = weighted_avg( Y_data )

    n = length(Y_data);
    
    if n > 0
    
    R = abs(Y_data'*ones(1,n)-ones(n,1)*Y_data);
    
    R_mean = mean(sum(R));
    
%     w = R_mean./sum(R);

    w = R_mean.^2./sum(R).^2;
    
    w_avg = (w*Y_data')/sum(w);
    
    w_var = sum(w.^2*((Y_data-w_avg).^2)')/n;
    
    w_std = sqrt(w_var);
    
    else
        
        w_avg = 0;

        w_var = 0;

        w_std = 0;
        
    end

end

