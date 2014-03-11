function [X_max n_Xmax] = calc_max_abs_value(X,n_start,n_stop)

for i = 1:length(n_start)
    if isnan(n_start(i)) == 0 && isnan(n_stop(i)) == 0
        X_sub = X(:,i);
        Xabs_sub = abs(X_sub);
        
        Xabs_max(i,1) = max(Xabs_sub(n_start(i):n_stop(i)));
        n_Xmax(i,1) = round(mean(find(Xabs_sub==Xabs_max(i,1))));
        
        if isnan(n_Xmax) == 0
            X_max(i,1) = X_sub(n_Xmax(i,1));
        else
            X_max(i,1) = nan;
        end
    else
        X_max(i,1) = nan;
        n_Xmax(i,1) = nan;
    end
end

