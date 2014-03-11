function [X_max n_Xmax] = calc_max_value(X,n_start,n_stop)

for i = 1:length(n_start)
    if isnan(n_start(i)) == 0 && isnan(n_stop(i)) == 0
        X_sub = X(:,i);
        X_max(i,1) = max(X_sub(n_start(i):n_stop(i)));
        n_Xmax(i,1) = round(mean(find(X_sub==X_max(i,1))));
    else
        X_max(i,1) = nan;
        n_Xmax(i,1) = nan;
    end
end

