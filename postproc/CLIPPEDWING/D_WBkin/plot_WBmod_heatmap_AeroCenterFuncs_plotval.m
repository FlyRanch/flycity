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

t_hist = [t_wb_AeroCenterFuncs_bins];
t_hist = t_hist(:);
binx = 0:1/(nx-1):1;

biny_min = -90;
biny_max = 90;
biny_pitch_min = -180;
biny_pitch_max = 180;
biny_dev_min = -45;
biny_dev_max = 180;

t_bins = [t_wb_AeroCenterFuncs_bins];           
binx = t_bins(:,1);

%% Left wing = intact wing: AeroCenterFuncIntact fourier series

% stroke
var = (AeroCenterFuncIntact_plot-AeroCenterFuncIntact_steady)*[strokeMOD_wb_L_AeroCenterFuncs_bins];

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
% ylabel('Intact wing','fontsize',10) 

% rotation
var =(AeroCenterFuncIntact_plot-AeroCenterFuncIntact_steady)*[pitchMOD_wb_L_AeroCenterFuncs_bins];

var = -var(:);
biny = biny_pitch_min: (biny_pitch_max - biny_pitch_min)/ny :biny_pitch_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,2)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_pitch_max -biny_pitch_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_pitch_max:(biny_pitch_max-biny_pitch_min)/2:-biny_pitch_min],'YTicklabel',-[-biny_pitch_max:(biny_pitch_max-biny_pitch_min)/2:-biny_pitch_min],'fontsize',8) 
set(gca,'YTick',[-biny_pitch_max:(biny_pitch_max-biny_pitch_min)/2:-biny_pitch_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Intact wing','fontsize',10) 

% deviation
var = (AeroCenterFuncIntact_plot-AeroCenterFuncIntact_steady)*[devMOD_wb_L_AeroCenterFuncs_bins];

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
ylabel('Intact wing','fontsize',10) 

% plot mean & std OR CI
subplot(3,3,1)
var = [strokeMOD_wb_L_AeroCenterFuncs_bins];
% plot_WBmeanstd
% plot_WBmeanCI
% plot_WBmeanCI_AeroCenterFuncIntact
plot_WBmeanCI_AeroCenterFuncIntact_NOcirc

subplot(3,3,2)
var = [pitchMOD_wb_L_AeroCenterFuncs_bins];
% plot_WBmeanstd
% plot_WBmeanCI
% plot_WBmeanCI_AeroCenterFuncIntact
plot_WBmeanCI_AeroCenterFuncIntact_NOcirc

subplot(3,3,3)
var = [devMOD_wb_L_AeroCenterFuncs_bins];
% plot_WBmeanstd
% plot_WBmeanCI
% plot_WBmeanCI_AeroCenterFuncIntact
plot_WBmeanCI_AeroCenterFuncIntact_NOcirc

% fourier series
subplot(3,3,1)
val=feval(strokeMOD_L_AeroCenterFuncs_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],(AeroCenterFuncIntact_plot-AeroCenterFuncIntact_steady)*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel('Intact wing','fontsize',10) 

subplot(3,3,2)
val=feval(pitchMOD_L_AeroCenterFuncs_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],(AeroCenterFuncIntact_plot-AeroCenterFuncIntact_steady)*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,3)
val=feval(devMOD_L_AeroCenterFuncs_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],(AeroCenterFuncIntact_plot-AeroCenterFuncIntact_steady)*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

%% Right wing = clipped wing: AeroCenterFuncClipped fourier series

% stroke
var = (AeroCenterFuncClipped_plot-AeroCenterFuncClipped_steady)*[strokeMOD_wb_R_AeroCenterFuncs_bins];

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
ylabel('Clipped wing','fontsize',10) 

% wingpitch
var = (AeroCenterFuncClipped_plot-AeroCenterFuncClipped_steady)*[pitchMOD_wb_R_AeroCenterFuncs_bins];

var = -var(:);
biny = biny_pitch_min: (biny_pitch_max - biny_pitch_min)/ny :biny_pitch_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,5)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_pitch_max -biny_pitch_min])
%     set(gca,'XTick',0:.5:1) 
    set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_pitch_max:(biny_pitch_max-biny_pitch_min)/2:-biny_pitch_min],'YTicklabel',-[-biny_pitch_max:(biny_pitch_max-biny_pitch_min)/2:-biny_pitch_min],'fontsize',8) 
set(gca,'YTick',[-biny_pitch_max:(biny_pitch_max-biny_pitch_min)/2:-biny_pitch_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Clipped wing','fontsize',10) 

% Deviation
var = (AeroCenterFuncClipped_plot-AeroCenterFuncClipped_steady)*[devMOD_wb_R_AeroCenterFuncs_bins];

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
ylabel('Clipped wing','fontsize',10) 

% plot mean & std OR CI
subplot(3,3,4)
var = [strokeMOD_wb_R_AeroCenterFuncs_bins];
% plot_WBmeanstd
% plot_WBmeanCI
% plot_WBmeanCI_AeroCenterFuncClipped
plot_WBmeanCI_AeroCenterFuncClipped_NOcirc

subplot(3,3,5)
var = [pitchMOD_wb_R_AeroCenterFuncs_bins];
% plot_WBmeanstd
% plot_WBmeanCI
% plot_WBmeanCI_AeroCenterFuncClipped
plot_WBmeanCI_AeroCenterFuncClipped_NOcirc

subplot(3,3,6)
var = [devMOD_wb_R_AeroCenterFuncs_bins];
% plot_WBmeanstd
% plot_WBmeanCI
% plot_WBmeanCI_AeroCenterFuncClipped
plot_WBmeanCI_AeroCenterFuncClipped_NOcirc


% fourier series
subplot(3,3,4)
val=feval(strokeMOD_R_AeroCenterFuncs_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],(AeroCenterFuncClipped_plot-AeroCenterFuncClipped_steady)*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel('Clipped wing','fontsize',10) 

subplot(3,3,5)
val=feval(pitchMOD_R_AeroCenterFuncs_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],(AeroCenterFuncClipped_plot-AeroCenterFuncClipped_steady)*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,6)
val=feval(devMOD_R_AeroCenterFuncs_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],(AeroCenterFuncClipped_plot-AeroCenterFuncClipped_steady)*val,'r')
legend off
xlabel([],'fontsize',10) 
ylabel([],'fontsize',10) 

%% Left - Right (intact - clipped) wingbeat kin:  AeroCenterFuncRatio

% stroke
var = (AeroCenterFuncRatio_plot-AeroCenterFuncRatio_steady)*[DstrokeMOD_wb_AeroCenterFuncs_bins];

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
ylabel('Clipped wing','fontsize',10) 

% Rotation
var = (AeroCenterFuncRatio_plot-AeroCenterFuncRatio_steady)*[DpitchMOD_wb_AeroCenterFuncs_bins];

var = -var(:);
biny = biny_pitch_min: (biny_pitch_max - biny_pitch_min)/ny :biny_pitch_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,8)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_pitch_max -biny_pitch_min])
    set(gca,'XTick',0:.5:1) 
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_pitch_max:(biny_pitch_max-biny_pitch_min)/2:-biny_pitch_min],'YTicklabel',-[-biny_pitch_max:(biny_pitch_max-biny_pitch_min)/2:-biny_pitch_min],'fontsize',8) 
set(gca,'YTick',[-biny_pitch_max:(biny_pitch_max-biny_pitch_min)/2:-biny_pitch_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Clipped wing','fontsize',10) 

% deviation
var = (AeroCenterFuncRatio_plot-AeroCenterFuncRatio_steady)*[DdevMOD_wb_AeroCenterFuncs_bins];

var = -var(:);
biny = biny_dev_min: (biny_dev_max - biny_dev_min)/ny :biny_dev_max;
var_hist = hist3([var,t_hist], {biny binx});
for i = 1:nx
    var_hist(:,i) = var_hist(:,i) / max(var_hist(:,i)); % normalize per time bin
end

subplot(3,3,9)
imagesc(binx,-biny,var_hist)
axis([0 1 -biny_dev_max -biny_dev_min])
    set(gca,'XTick',0:.5:1) 
%     set(gca,'XTick',0:.5:1,'XTickLabel',[]) 
% set(gca,'YTick',[-biny_dev_max:(biny_dev_max-biny_dev_min)/2:-biny_dev_min],'YTicklabel',-[-biny_dev_max:(biny_dev_max-biny_dev_min)/2:-biny_dev_min],'fontsize',8) 
set(gca,'YTick',[-biny_dev_max:(biny_dev_max-biny_dev_min)/2:-biny_dev_min],'fontsize',8) 
% xlabel('time','fontsize',10) 
ylabel('Clipped wing','fontsize',10) 

%% plot mean & std OR CI
subplot(3,3,7)
var = [DstrokeMOD_wb_AeroCenterFuncs_bins];
% plot_WBmeanstd
% plot_WBmeanCI
% plot_WBmeanCI_AeroCenterFuncRatio
plot_WBmeanCI_AeroCenterFuncRatio_NOcirc

subplot(3,3,8)
var = [DpitchMOD_wb_AeroCenterFuncs_bins];
% plot_WBmeanstd
% plot_WBmeanCI
% plot_WBmeanCI_AeroCenterFuncRatio
plot_WBmeanCI_AeroCenterFuncRatio_NOcirc

subplot(3,3,9)
var = [DdevMOD_wb_AeroCenterFuncs_bins];
% plot_WBmeanstd
% plot_WBmeanCI
% plot_WBmeanCI_AeroCenterFuncRatio
plot_WBmeanCI_AeroCenterFuncRatio_NOcirc

%% fourier series
subplot(3,3,7)
val=feval(DstrokeMOD_AeroCenterFuncs_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],(AeroCenterFuncRatio_plot-AeroCenterFuncRatio_steady)*val,'r')
legend off
xlabel('time','fontsize',10) 
ylabel('Intact - Clipped','fontsize',10) 

subplot(3,3,8)
val=feval(DpitchMOD_AeroCenterFuncs_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],(AeroCenterFuncRatio_plot-AeroCenterFuncRatio_steady)*val,'r')
legend off
xlabel('time','fontsize',10) 
ylabel([],'fontsize',10) 

subplot(3,3,9)
val=feval(DdevMOD_AeroCenterFuncs_fourier_fit_binmean,[0:.01:1]);
plot([0:.01:1],(AeroCenterFuncRatio_plot-AeroCenterFuncRatio_steady)*val,'r')
legend off
xlabel('time','fontsize',10) 
ylabel([],'fontsize',10) 


