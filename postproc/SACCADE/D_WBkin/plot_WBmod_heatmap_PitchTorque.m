figure
subplot(3,1,1)
title('Pitch Torque wing mod','fontsize',12)
hold on
subplot(3,1,2)
% title('wing pitch')
hold on
subplot(3,1,3)
% title('stroke deviation')
hold on

colormap(cmap)
color_band = 'b';

% bins
nx=n;
ny=100;

t_hist = [t_wb_PitchTorque_bins t_wb_PitchTorque_bins];
t_hist = t_hist(:);
binx = 0:1/(nx-1):1;

biny_min = -45;
biny_max = 45;

%% stroke
% biny_min = -30;
% biny_max = 30;
var = [strokeMOD_wb_L_PitchTorque_bins strokeMOD_wb_R_PitchTorque_bins];

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
% ylabel('Left wing','fontsize',10) 
% title('wing stroke')

%% wingpitchMOD
% biny_min = -30;
% biny_max = 30;
var = [pitchMOD_wb_L_PitchTorque_bins pitchMOD_wb_R_PitchTorque_bins];

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
% title('wing pitch')

%% devMOD
% biny_min = -30;
% biny_max = 30;
var = [devMOD_wb_L_PitchTorque_bins devMOD_wb_R_PitchTorque_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,1,3)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
    set(gca,'XTick',0:.5:1) 
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('Left wing','fontsize',10) 
% title('stroke deviation')

%% plot mean & std OR CI
binx = t_wb_PitchTorque_bins(:,1);
t_bins = [t_wb_PitchTorque_bins t_wb_PitchTorque_bins];           

subplot(3,1,1)
var = [strokeMOD_wb_L_PitchTorque_bins strokeMOD_wb_R_PitchTorque_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,1,2)
var = [pitchMOD_wb_L_PitchTorque_bins pitchMOD_wb_R_PitchTorque_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,1,3)
var = [devMOD_wb_L_PitchTorque_bins devMOD_wb_R_PitchTorque_bins];
% plot_WBmeanstd
plot_WBmeanCI

%% plot polinomials
% subplot(3,1,1)
% plot_WBfitting_singlevar_updowncut(strokeMOD_PitchTorque_fit_binmean_periodic,'r','g');
% % plot_WBfitting_singlevar_updowncut(strokeMOD_PitchTorque_fit_binmean,'g');
% 
% subplot(3,1,2)
% plot_WBfitting_singlevar_updowncut(pitchMOD_PitchTorque_fit_binmean_periodic,'r','g');
% % plot_WBfitting_singlevar_updowncut(pitchMOD_PitchTorque_fit_binmean,'g');
% 
% subplot(3,1,3)
% plot_WBfitting_singlevar_updowncut(devMOD_PitchTorque_fit_binmean_periodic,'r','g');
% % plot_WBfitting_singlevar_updowncut(devMOD_PitchTorque_fit_binmean,'g');


%% polynomials
subplot(3,1,1)
calc_val_fourier_series_4thN8th_order(binx,strokeMOD_PitchTorque_fourier_coeffs_binmean,plot_fits);
% plot(strokeMOD_PitchTorque_fourier_fit_binmean);
legend off
xlabel([],'fontsize',10) 
ylabel('wing stroke','fontsize',10) 

subplot(3,1,2)
calc_val_fourier_series_4thN8th_order(binx,pitchMOD_PitchTorque_fourier_coeffs_binmean,plot_fits);
% plot(pitchMOD_PitchTorque_fourier_fit_binmean);
legend off
xlabel([],'fontsize',10) 
ylabel('wing pitch','fontsize',10) 

subplot(3,1,3)
calc_val_fourier_series_4thN8th_order(binx,devMOD_PitchTorque_fourier_coeffs_binmean,plot_fits);
% plot(devMOD_PitchTorque_fourier_fit_binmean);
legend off
xlabel('normalized time','fontsize',10) 
ylabel('wing deviation','fontsize',10) 




