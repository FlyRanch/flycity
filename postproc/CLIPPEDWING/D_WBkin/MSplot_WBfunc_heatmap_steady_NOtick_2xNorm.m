figure
subplot(1,3,1)
% title('steady wingbeat','fontsize',12)
hold on
subplot(1,3,2)
% title('wing pitch')
hold on
subplot(1,3,3)
% title('stroke deviation')
hold on

colormap(cmap)
color_band = [.5 .5 .5];

% bins
nx=n;
ny=200;

nx=200;
ny=300;

t_hist = [t_wb_steady_bins t_wb_steady_bins];
t_hist = t_hist(:);
binx = 0:1/(nx-1):1;

%% stroke
biny_min = -90;
biny_max = 90;
var = [stroke_wb_L_steady_bins stroke_wb_L_steady_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(1,3,1)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'XTick',[]) 
set(gca,'YTick',[]) 
% xlabel('time','fontsize',10) 
% ylabel('Left wing','fontsize',10) 

%% wingpitch
biny_min = -180;
biny_max = 0;
var = [pitch_wb_L_steady_bins pitch_wb_L_steady_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(1,3,2)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'XTick',[]) 
set(gca,'YTick',[]) 
% xlabel('time','fontsize',10) 
% ylabel('Left wing','fontsize',10) 

%% dev
biny_min = -30;
biny_max = 30;
var = [dev_wb_L_steady_bins dev_wb_L_steady_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(1,3,3)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'XTick',[]) 
set(gca,'YTick',[]) 
% xlabel('time','fontsize',10) 
% ylabel('Left wing','fontsize',10) 

%% mean & std            
% binx = t_wb_steady_bins(:,1);
% t_bins = [t_wb_steady_bins t_wb_steady_bins];           
% 
% 
%             subplot(1,3,1)
%             var = [stroke_wb_L_steady_bins stroke_wb_R_steady_bins];
% %             plot_WBmeanstd
%             plot_WBmeanCI
%             
%             subplot(1,3,2)
%             var = [pitch_wb_L_steady_bins pitch_wb_R_steady_bins];
% %             plot_WBmeanstd
%             plot_WBmeanCI
% 
%             subplot(1,3,3)
%             var = [dev_wb_L_steady_bins dev_wb_R_steady_bins];
% %             plot_WBmeanstd
%             plot_WBmeanCI

%% polynomials
% subplot(1,3,1)
% plot_WBfitting_singlevar_updowncut(stroke_steady_fit_binmean_periodic,'r','g');
% % plot_WBfitting_singlevar_updowncut(stroke_steady_fit_binmean,'g','g');
% 
% subplot(1,3,2)
% plot_WBfitting_singlevar_updowncut(pitch_steady_fit_binmean_periodic,'r','g');
% % plot_WBfitting_singlevar_updowncut(pitch_steady_fit_binmean,'g','g');
% 
% subplot(1,3,3)
% plot_WBfitting_singlevar_updowncut(dev_steady_fit_binmean_periodic,'r','g');
% % plot_WBfitting_singlevar_updowncut(dev_steady_fit_binmean,'g','g');


%% fourier series
% subplot(1,3,1)
% calc_val_fourier_series_4thN8th_order(binx,stroke_steady_fourier_coeffs_binmean,plot_fits);
% % plot(stroke_steady_fourier_fit_binmean);
% legend off
% xlabel([],'fontsize',10) 
% ylabel('wing stroke','fontsize',10) 
% 
% subplot(1,3,2)
% calc_val_fourier_series_4thN8th_order(binx,pitch_steady_fourier_coeffs_binmean,plot_fits);
% % plot(pitch_steady_fourier_fit_binmean);
% legend off
% xlabel([],'fontsize',10) 
% ylabel('wing pitch','fontsize',10) 
% 
% subplot(1,3,3)
% calc_val_fourier_series_4thN8th_order(binx,dev_steady_fourier_coeffs_binmean,plot_fits);
% % plot(dev_steady_fourier_fit_binmean);
% legend off
% xlabel('normalized time','fontsize',10) 
% ylabel('wing deviation','fontsize',10) 
% 
% 



