for i = 1:length(binx)

    var_now = var(t_bins == binx(i));
    [mu ul ll] = circ_mean_deg_nonan(var_now);
    [std std0] = circ_std_deg_nonan(var_now);
    var_mean(i,1) = mu;
    var_ul(i,1) = ul;
    var_ll(i,1) = ll;
    var_meanCI(i,:) = [mu ul ll];
    
    var_meanstd(i,:) = [mu mu+std mu-std];
%     var_ul(i,1) = mu+std;
%     var_ll(i,1) = mu-std;
end

ciplot(var_ll,var_ul,binx,color_band)
% alpha(.25)
hold on
plot(binx,var_mean,'k')
% ylim([y_min,y_max])
% set(gca,'XTick',[0;0.5;1],'fontsize',8)  
% set(gca,'YTick',[y_min;0;y_max],'fontsize',8)  
% grid on
