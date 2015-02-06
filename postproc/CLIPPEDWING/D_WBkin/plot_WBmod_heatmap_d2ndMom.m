figure
subplot(3,3,1)
title('wing stroke')
hold on
subplot(3,3,2)
title('wing pitch')
hold on
subplot(3,3,3)
title('stroke deviation')
hold on
subplot(3,3,4)
hold on
subplot(3,3,5)
hold on
subplot(3,3,6)
hold on
subplot(3,3,7)
hold on
subplot(3,3,8)
hold on
subplot(3,3,9)
hold on

colormap(cmap)
color_band = 'b';

% bins
nx=n;
ny=100;

t_hist = [t_wb_d2ndMom_bins];
t_hist = t_hist(:);
binx = 0:1/(nx-1):1;

biny_min = -45;
biny_max = 45;

%% stroke

% in rotRaccel Left wing
% biny_min = -30;
% biny_max = 30;
var = [strokeMOD_wb_L_d2ndMom_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,1)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
% ylabel('Left wing','fontsize',10) 

% in rotRaccel Right wing
% biny_min = -30;
% biny_max = 30;
var = [strokeMOD_wb_R_d2ndMom_bins];

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
ylabel('Right wing','fontsize',10) 

% for-rear
% biny_min = -30;
% biny_max = 30;
var = [DstrokeMOD_wb_d2ndMom_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,7)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
    set(gca,'XTick',0:.5:1) 
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Right wing','fontsize',10) 


%% wingpitch

% in rotRaccel Left wing
% biny_min = -30;
% biny_max = 30;
var = [pitchMOD_wb_L_d2ndMom_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,2)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Left wing','fontsize',10) 

% in rotRaccel Right wing
% biny_min = -30;
% biny_max = 30;
var = [pitchMOD_wb_R_d2ndMom_bins];

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
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Right wing','fontsize',10) 

% for-rear
% biny_min = -30;
% biny_max = 30;
var = [DpitchMOD_wb_d2ndMom_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,8)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
    set(gca,'XTick',0:.5:1) 
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Right wing','fontsize',10) 


%% dev

% in rotRaccel Left wing
% biny_min = -30;
% biny_max = 30;
var = [devMOD_wb_L_d2ndMom_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,3)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Left wing','fontsize',10) 

% in rotRaccel Right wing
% biny_min = -30;
% biny_max = 30;
var = [devMOD_wb_R_d2ndMom_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,6)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Right wing','fontsize',10) 

% for-rear
% biny_min = -30;
% biny_max = 30;
var = [DdevMOD_wb_d2ndMom_bins];

var = -var(:);
biny = biny_min: (biny_max - biny_min)/ny :biny_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,9)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_max -biny_min])
    set(gca,'XTick',0:.5:1) 
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'YTicklabel',-[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
set(gca,'YTick',[-biny_max:(biny_max-biny_min)/2:-biny_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Right wing','fontsize',10) 


%% plot mean & std OR CI
t_bins = [t_wb_d2ndMom_bins];           
binx = t_bins(:,1);

subplot(3,3,1)
var = [strokeMOD_wb_L_d2ndMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,2)
var = [pitchMOD_wb_L_d2ndMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,3)
var = [devMOD_wb_L_d2ndMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,4)
var = [strokeMOD_wb_R_d2ndMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,5)
var = [pitchMOD_wb_R_d2ndMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,6)
var = [devMOD_wb_R_d2ndMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,7)
var = [DstrokeMOD_wb_d2ndMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,8)
var = [DpitchMOD_wb_d2ndMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,9)
var = [DdevMOD_wb_d2ndMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

%% plot polinomials
% subplot(3,3,1)
% plot_WBfitting_singlevar_updowncut(strokeMOD_L_d2ndMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,2)
% plot_WBfitting_singlevar_updowncut(pitchMOD_L_d2ndMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,3)
% plot_WBfitting_singlevar_updowncut(devMOD_L_d2ndMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,4)
% plot_WBfitting_singlevar_updowncut(strokeMOD_R_d2ndMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,5)
% plot_WBfitting_singlevar_updowncut(pitchMOD_R_d2ndMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,6)
% plot_WBfitting_singlevar_updowncut(devMOD_R_d2ndMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,7)
% plot_WBfitting_singlevar_updowncut(DstrokeMOD_d2ndMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,8)
% plot_WBfitting_singlevar_updowncut(DpitchMOD_d2ndMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,9)
% plot_WBfitting_singlevar_updowncut(DdevMOD_d2ndMom_fit_binmean_periodic,'r','g');

%% fourier series

subplot(3,3,1)
plot(strokeMOD_L_d2ndMom_fourier_fit_binmean);
legend off
xlabel([],'fontsize',10) 
ylabel('Left wing','fontsize',10) 

subplot(3,3,2)
plot(pitchMOD_L_d2ndMom_fourier_fit_binmean);
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,3)
plot(devMOD_L_d2ndMom_fourier_fit_binmean);
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,4)
plot(strokeMOD_R_d2ndMom_fourier_fit_binmean);
legend off
xlabel([],'fontsize',10) 
ylabel('Right wing','fontsize',10) 

subplot(3,3,5)
plot(pitchMOD_R_d2ndMom_fourier_fit_binmean);
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,6)
plot(devMOD_R_d2ndMom_fourier_fit_binmean);
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,7)
plot(DstrokeMOD_d2ndMom_fourier_fit_binmean);
legend off
xlabel('time','fontsize',10) 
ylabel('Left - Right','fontsize',10) 

subplot(3,3,8)
plot(DpitchMOD_d2ndMom_fourier_fit_binmean);
legend off
xlabel('time','fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,9)
plot(DdevMOD_d2ndMom_fourier_fit_binmean);
legend off
xlabel('time','fontsize',10) 
ylabel([],'fontsize',10) 





