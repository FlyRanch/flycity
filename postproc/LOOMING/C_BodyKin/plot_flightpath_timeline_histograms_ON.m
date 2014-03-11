%% plot flightpath timeline histograms ON

t_hist = repmat(t,136,1);

[N_histx,x_hist] = hist(t_hist,100);
[N_histy, y_hist] = hist(var,360);

var_hist = hist3([t_hist, var], [100,360]);

    figure()
    image(y_hist,x_hist,var_hist)
