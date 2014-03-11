

% figure
subplot(3,3,4)
% title('all wingbeats','fontsize',12)
hold on
subplot(3,3,5)
% title('wing pitch')
hold on
subplot(3,3,6)
% title('stroke deviation')
hold on

colormap(cmap)
color_band = [.5 .5 .5];

% bins
nx=n;
ny=100;

t_hist = [t_wb_bins t_wb_bins];
t_hist = t_hist(:);
binx = 0:1/(nx-1):1;

%% stroke
biny_min = -90;
biny_max = 90;
var = [stroke_wb_L_bins stroke_wb_L_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,4)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('Left wing','fontsize',10) 

%% wingpitch
biny_min = -180;
biny_max = 0;
var = [pitch_wb_L_bins pitch_wb_L_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,5)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('Left wing','fontsize',10) 

%% dev
biny_min = -30;
biny_max = 30;
var = [dev_wb_L_bins dev_wb_L_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,6)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
    set(gca,'XTick',0:.5:1) 
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('Left wing','fontsize',10) 

%% mean & std            
% binx = t_wb_bins(:,1);
% t_bins = [t_wb_bins t_wb_bins];           
% 
% 
%             subplot(3,3,4)
%             var = [stroke_wb_L_bins stroke_wb_R_bins];
% %             plot_WBmeanstd
%             plot_WBmeanCI
%             
%             subplot(3,3,5)
%             var = [pitch_wb_L_bins pitch_wb_R_bins];
% %             plot_WBmeanstd
%             plot_WBmeanCI
% 
%             subplot(3,3,6)
%             var = [dev_wb_L_bins dev_wb_R_bins];
% %             plot_WBmeanstd
%             plot_WBmeanCI

%% polynomials
% subplot(3,3,4)
% plot_WBfitting_singlevar_updowncut(stroke_fit_binmean_periodic,'r','g');
% % plot_WBfitting_singlevar_updowncut(stroke_fit_binmean,'g','g');
% 
% subplot(3,3,5)
% plot_WBfitting_singlevar_updowncut(pitch_fit_binmean_periodic,'r','g');
% % plot_WBfitting_singlevar_updowncut(pitch_fit_binmean,'g','g');
% 
% subplot(3,3,6)
% plot_WBfitting_singlevar_updowncut(dev_fit_binmean_periodic,'r','g');
% % plot_WBfitting_singlevar_updowncut(dev_fit_binmean,'g','g');


%% fourier series
% subplot(3,3,4)
% calc_val_fourier_series_4thN8th_order(binx,stroke_fourier_coeffs_binmean,plot_fits);
% % plot(stroke_fourier_fit_binmean);
% legend off
% xlabel([],'fontsize',10) 
% ylabel('wing stroke','fontsize',10) 
% 
% subplot(3,3,5)
% calc_val_fourier_series_4thN8th_order(binx,pitch_fourier_coeffs_binmean,plot_fits);
% % plot(pitch_fourier_fit_binmean);
% legend off
% xlabel([],'fontsize',10) 
% ylabel('wing pitch','fontsize',10) 
% 
% subplot(3,3,6)
% calc_val_fourier_series_4thN8th_order(binx,dev_fourier_coeffs_binmean,plot_fits);
% % plot(dev_fourier_fit_binmean);
% legend off
% xlabel('normalized time','fontsize',10) 
% ylabel('wing deviation','fontsize',10) 





