% bins
nx=n;
ny=2*100;
binx = 0:1/(nx-1):1;

% time
t_hist_steady = [t_wb_steady_bins t_wb_steady_bins];
t_hist_steady = t_hist_steady(:);
t_hist_NONsteady = [t_wb_NONsteady_bins t_wb_NONsteady_bins];
t_hist_NONsteady = t_hist_NONsteady(:);

%% stroke
biny_min = -90;
biny_max = 90;
biny = biny_min: (biny_max - biny_min)/ny :biny_max;

var_steady = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
var_steady = -var_steady(:);
var_hist_steady = hist3([var_steady,t_hist_steady], {biny binx});

var_NONsteady = [stroke_wb_L_NONsteady_bins stroke_wb_R_NONsteady_bins];
var_NONsteady = -var_NONsteady(:);
var_hist_NONsteady = hist3([var_NONsteady,t_hist_NONsteady], {biny binx});

for i = 1:nx
    var_hist_steady(:,i) = var_hist_steady(:,i) / max(var_hist_steady(:,i)); % normalize per time bin
    var_hist_NONsteady(:,i) = var_hist_NONsteady(:,i) / max(var_hist_NONsteady(:,i)); % normalize per time bin
end

% inverse
var_hist_steady = 1 - var_hist_steady;
var_hist_NONsteady = 1 - var_hist_NONsteady;

clear var_hist
var_hist(:,:) = var_hist_NONsteady;
var_hist(:,:,2) = var_hist_steady;
var_hist(:,:,3) = var_hist_steady;


figure
imagesc(binx,biny,var_hist)
axis([0 1 biny_min biny_max])
    set(gca,'XTick',[],'XTickLabel',[]) 
set(gca,'YTick',[],'fontsize',8) 
% axis off

imwrite(var_hist,['MSplot_WBsteadyNmod_heatmap_stroke.tif'])
saveas(gca,['MSplot_WBsteadyNmod_heatmap_stroke.fig'])
saveas(gca,['MSplot_WBsteadyNmod_heatmap_stroke.png'])
plot2svg(['MSplot_WBsteadyNmod_heatmap_stroke.svg'])    


%% pitch
biny_min = -180;
biny_max = 0;
biny = biny_min: (biny_max - biny_min)/ny :biny_max;

var_steady = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
var_steady = -var_steady(:);
var_hist_steady = hist3([var_steady,t_hist_steady], {biny binx});

var_NONsteady = [pitch_wb_L_NONsteady_bins pitch_wb_R_NONsteady_bins];
var_NONsteady = -var_NONsteady(:);
var_hist_NONsteady = hist3([var_NONsteady,t_hist_NONsteady], {biny binx});

for i = 1:nx
    var_hist_steady(:,i) = var_hist_steady(:,i) / max(var_hist_steady(:,i)); % normalize per time bin
    var_hist_NONsteady(:,i) = var_hist_NONsteady(:,i) / max(var_hist_NONsteady(:,i)); % normalize per time bin
end

% inverse
var_hist_steady = 1 - var_hist_steady;
var_hist_NONsteady = 1 - var_hist_NONsteady;

clear var_hist
var_hist(:,:) = var_hist_NONsteady;
var_hist(:,:,2) = var_hist_steady;
var_hist(:,:,3) = var_hist_steady;


figure
imagesc(binx,biny,var_hist)
axis([0 1 biny_min biny_max])
    set(gca,'XTick',[],'XTickLabel',[]) 
set(gca,'YTick',[],'fontsize',8) 
% axis off

imwrite(var_hist,['MSplot_WBsteadyNmod_heatmap_pitch.tif'])
saveas(gca,['MSplot_WBsteadyNmod_heatmap_pitch.fig'])
saveas(gca,['MSplot_WBsteadyNmod_heatmap_pitch.png'])
plot2svg(['MSplot_WBsteadyNmod_heatmap_pitch.svg'])    


%% dev
biny_min = -30;
biny_max = 30;
biny = biny_min: (biny_max - biny_min)/ny :biny_max;

var_steady = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
var_steady = -var_steady(:);
var_hist_steady = hist3([var_steady,t_hist_steady], {biny binx});

var_NONsteady = [dev_wb_L_NONsteady_bins dev_wb_R_NONsteady_bins];
var_NONsteady = -var_NONsteady(:);
var_hist_NONsteady = hist3([var_NONsteady,t_hist_NONsteady], {biny binx});

for i = 1:nx
    var_hist_steady(:,i) = var_hist_steady(:,i) / max(var_hist_steady(:,i)); % normalize per time bin
    var_hist_NONsteady(:,i) = var_hist_NONsteady(:,i) / max(var_hist_NONsteady(:,i)); % normalize per time bin
end

% inverse
var_hist_steady = 1 - var_hist_steady;
var_hist_NONsteady = 1 - var_hist_NONsteady;

clear var_hist
var_hist(:,:) = var_hist_NONsteady;
var_hist(:,:,2) = var_hist_steady;
var_hist(:,:,3) = var_hist_steady;


figure
imagesc(binx,biny,var_hist)
axis([0 1 biny_min biny_max])
    set(gca,'XTick',[],'XTickLabel',[]) 
set(gca,'YTick',[],'fontsize',8) 
% axis off

imwrite(var_hist,['MSplot_WBsteadyNmod_heatmap_dev.tif'])
saveas(gca,['MSplot_WBsteadyNmod_heatmap_dev.fig'])
saveas(gca,['MSplot_WBsteadyNmod_heatmap_dev.png'])
plot2svg(['MSplot_WBsteadyNmod_heatmap_dev.svg'])    



