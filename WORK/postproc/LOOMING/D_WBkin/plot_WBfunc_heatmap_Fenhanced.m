figure
subplot(3,1,1)
title('wing stroke')
hold on
subplot(3,1,2)
title('wing pitch')
hold on
subplot(3,1,3)
title('stroke deviation')
hold on

colormap(cmap)

% bins
nx=n;
ny=100;

t_hist = [t_wb_Fenhance_bins t_wb_Fenhance_bins];
t_hist = t_hist(:);
binx = 0:1/(nx-1):1;

% stroke
biny_min = -90;
biny_max = 90;
var = [stroke_wb_L_Fenhance_bins stroke_wb_L_Fenhance_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,1,1)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Left wing','fontsize',10) 
title('wing stroke')

% wingpitch
biny_min = -180;
biny_max = 0;
var = [pitch_wb_L_Fenhance_bins pitch_wb_L_Fenhance_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,1,2)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('Left wing','fontsize',10) 
title('wing pitch')

% dev
biny_min = -30;
biny_max = 30;
var = [dev_wb_L_Fenhance_bins dev_wb_L_Fenhance_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,1,3)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('Left wing','fontsize',10) 
title('stroke deviation')

% plot_WBfunc_csaps_LRdiff
% plot_WBfunc_circmeanCI
plot_WBfunc_Fenhanced_circmeanstd


subplot(3,1,1)
plot_WBfitting_singlevar_updowncut(stroke_Fenhance_fit_binmean_periodic,'r','g');
% plot_WBfitting_singlevar_updowncut(stroke_Fenhance_fit_binmean,'g');

subplot(3,1,2)
plot_WBfitting_singlevar_updowncut(pitch_Fenhance_fit_binmean_periodic,'r','g');
% plot_WBfitting_singlevar_updowncut(pitch_Fenhance_fit_binmean,'g');

subplot(3,1,3)
plot_WBfitting_singlevar_updowncut(dev_Fenhance_fit_binmean_periodic,'r','g');
% plot_WBfitting_singlevar_updowncut(dev_Fenhance_fit_binmean,'g');





