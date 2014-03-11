function X_mean = calc_mean_value(X,n_start,n_stop)

for i = 1:length(n_start)
    if isnan(n_start(i)) == 0 && isnan(n_stop(i)) == 0
        X_mean(i,1) = mean(X(n_start(i):n_stop(i),i));
    else
        X_mean(i,1) = nan;
    end
end

