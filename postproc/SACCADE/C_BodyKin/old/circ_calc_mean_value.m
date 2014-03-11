function [X_mean X_ul X_ll] = circ_calc_mean_value(X,n_start,n_stop)

for i = 1:length(n_start)
    if isnan(n_start(i)) == 0 && isnan(n_stop(i)) == 0
        [mu ul ll] = circ_mean_deg_nonan(X(n_start(i):n_stop(i),i));
        X_mean(i,1) = mu;
        X_ul(i,1) = ul;
        X_ll(i,1) = ll;
    else
        X_mean(i,1) = nan;
        X_ul(i,1) = nan;
        X_ll(i,1) = nan;
    end
end

