for i = 1:length(binx)

    var_now = var(t_bins == binx(i));
    mu = nanmean(var_now);
    std = nanstd(var_now);
    CI = 1.96 * std / sqrt(length(var_now));
    ul = mu + CI;
    ll = mu - CI;

    var_mean(i,1) = mu;
    var_ul(i,1) = ul;
    var_ll(i,1) = ll;
    var_meanCI(i,:) = [mu ul ll];
    
    var_meanstd(i,:) = [mu mu+std mu-std];
%     var_ul(i,1) = mu+std;
%     var_ll(i,1) = mu-std;
end

% % var_mean = radtodeg(unwrap(degtorad(var_mean)));
% % var_ll = radtodeg(unwrap(degtorad(var_ll)));
% % var_ul = radtodeg(unwrap(degtorad(var_ul)));
% 
% var_ll = naninterp(var_ll);
% var_ul = naninterp(var_ul);

ciplot((AeroCenterFuncRatio_plot-AeroCenterFuncRatio_steady)*var_ll,(AeroCenterFuncRatio_plot-AeroCenterFuncRatio_steady)*var_ul,binx,color_band)
% alpha(.25)
hold on
plot(binx,(AeroCenterFuncRatio_plot-AeroCenterFuncRatio_steady)*var_mean,'k')
% ylim([y_min,y_max])
% set(gca,'XTick',[0;0.5;1],'fontsize',8)  
% set(gca,'YTick',[y_min;0;y_max],'fontsize',8)  
% grid on
