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

t_hist = [t_wb_d1stMom_bins];
t_hist = t_hist(:);
binx = 0:1/(nx-1):1;

biny_min = -90;
biny_max = 90;

%% stroke

% in rotRaccel Left wing
% biny_min = -30;
% biny_max = 30;
var = FirstMomentRatio_plot*[strokeMOD_wb_L_d1stMom_bins];

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
var = FirstMomentRatio_plot*[strokeMOD_wb_R_d1stMom_bins];

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
var = FirstMomentRatio_plot*[DstrokeMOD_wb_d1stMom_bins];

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
var = FirstMomentRatio_plot*[pitchMOD_wb_L_d1stMom_bins];

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
var = FirstMomentRatio_plot*[pitchMOD_wb_R_d1stMom_bins];

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
var = FirstMomentRatio_plot*[DpitchMOD_wb_d1stMom_bins];

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
var = FirstMomentRatio_plot*[devMOD_wb_L_d1stMom_bins];

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
var = FirstMomentRatio_plot*[devMOD_wb_R_d1stMom_bins];

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
var = FirstMomentRatio_plot*[DdevMOD_wb_d1stMom_bins];

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
t_bins = [t_wb_d1stMom_bins];           
binx = t_bins(:,1);

subplot(3,3,1)
var = FirstMomentRatio_plot*[strokeMOD_wb_L_d1stMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,2)
var = FirstMomentRatio_plot*[pitchMOD_wb_L_d1stMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,3)
var = FirstMomentRatio_plot*[devMOD_wb_L_d1stMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,4)
var = FirstMomentRatio_plot*[strokeMOD_wb_R_d1stMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,5)
var = FirstMomentRatio_plot*[pitchMOD_wb_R_d1stMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,6)
var = FirstMomentRatio_plot*[devMOD_wb_R_d1stMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,7)
var = FirstMomentRatio_plot*[DstrokeMOD_wb_d1stMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,8)
var = FirstMomentRatio_plot*[DpitchMOD_wb_d1stMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

subplot(3,3,9)
var = FirstMomentRatio_plot*[DdevMOD_wb_d1stMom_bins];
% plot_WBmeanstd
plot_WBmeanCI

%% plot polinomials
% subplot(3,3,1)
% plot_WBfitting_singlevar_updowncut(strokeMOD_L_d1stMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,2)
% plot_WBfitting_singlevar_updowncut(pitchMOD_L_d1stMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,3)
% plot_WBfitting_singlevar_updowncut(devMOD_L_d1stMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,4)
% plot_WBfitting_singlevar_updowncut(strokeMOD_R_d1stMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,5)
% plot_WBfitting_singlevar_updowncut(pitchMOD_R_d1stMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,6)
% plot_WBfitting_singlevar_updowncut(devMOD_R_d1stMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,7)
% plot_WBfitting_singlevar_updowncut(DstrokeMOD_d1stMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,8)
% plot_WBfitting_singlevar_updowncut(DpitchMOD_d1stMom_fit_binmean_periodic,'r','g');
% 
% subplot(3,3,9)
% plot_WBfitting_singlevar_updowncut(DdevMOD_d1stMom_fit_binmean_periodic,'r','g');

%% fourier series

subplot(3,3,1)
val=feval(strokeMOD_L_d1stMom_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],FirstMomentRatio_plot*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel('Left wing','fontsize',10) 

subplot(3,3,2)
val=feval(pitchMOD_L_d1stMom_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],FirstMomentRatio_plot*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,3)
val=feval(devMOD_L_d1stMom_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],FirstMomentRatio_plot*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,4)
val=feval(strokeMOD_R_d1stMom_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],FirstMomentRatio_plot*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel('Right wing','fontsize',10) 

subplot(3,3,5)
val=feval(pitchMOD_R_d1stMom_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],FirstMomentRatio_plot*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,6)
val=feval(devMOD_R_d1stMom_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],FirstMomentRatio_plot*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,7)
val=feval(DstrokeMOD_d1stMom_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],FirstMomentRatio_plot*val,'r')
legend off
xlabel('time','fontsize',10) 
ylabel('Left - Right','fontsize',10) 

subplot(3,3,8)
val=feval(DpitchMOD_d1stMom_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],FirstMomentRatio_plot*val,'r')
legend off
xlabel('time','fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,9)
val=feval(DdevMOD_d1stMom_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],FirstMomentRatio_plot*val,'r')
legend off
xlabel('time','fontsize',10) 
ylabel([],'fontsize',10) 





