function [X_min n_Xmin] = calc_min_value(X,n_start,n_stop)

for i = 1:length(n_start)
    if isnan(n_start(i)) == 0 && isnan(n_stop(i)) == 0
        X_sub = X(:,i);
        X_min(i,1) = min(X_sub(n_start(i):n_stop(i)));
        n_Xmin(i,1) = round(mean(find(X_sub==X_min(i,1))));
    else
        X_min(i,1) = nan;
        n_Xmin(i,1) = nan;
    end
end

