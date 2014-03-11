%% heatmap hist log

yx_hist = hist3([y_hist,x_hist], {biny binx});
yx_hist_log = log(yx_hist);
imagesc(binx,-biny,yx_hist_log)
