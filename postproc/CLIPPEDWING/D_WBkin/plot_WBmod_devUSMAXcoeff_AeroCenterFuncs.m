%% dev max

figure
subplot(3,2,1)
title('devmax@upstroke')
hold on
subplot(3,2,2)
title('devmax@up mod coeff')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

% max up dev bin location
start = 50;
stop = 150;
var_loc = dev_wb_steady_bins_meanCIstd(start:stop,1);
n_loc = start+find(var_loc == max(var_loc));

%% left wing = intact
subplot(3,2,1)
plot((AeroCenterFuncIntact),dev_wb_L_bins(n_loc,:),'.g')
plot((AeroCenterFuncIntact_AeroCenterFuncs),dev_wb_L_AeroCenterFuncs_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_L_AeroCenterFuncs_bins_meanCIstd(n_loc,1);
x0 = AeroCenterFuncIntact_steady;
x1 = floor(10*min(AeroCenterFuncIntact_AeroCenterFuncs(:)))/10;
y0 = dev_wb_steady_bins_meanCIstd(n_loc);
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('AeroCenterFuncIntact','fontsize',10)
ylabel('Intact wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([x1 x0 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% coeff
subplot(3,2,2)
plot((AeroCenterFuncIntact_AeroCenterFuncs),devMOD_wb_L_AeroCenterFuncs_bins(n_loc,:),'.b')

RC = devMOD_wb_L_AeroCenterFuncs_bins_meanCIstd(n_loc,1);
plot([x0 x1],[RC RC],'-r','linewidth',2)
xlabel('AeroCenterFuncIntact','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -180;
ang_max = 540;
axis([x1 x0 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% right wing = clipped
subplot(3,2,3)
plot((AeroCenterFuncClipped),dev_wb_R_bins(n_loc,:),'.g')
plot(AeroCenterFuncClipped_AeroCenterFuncs,dev_wb_R_AeroCenterFuncs_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_R_AeroCenterFuncs_bins_meanCIstd(n_loc,1);
x0 = AeroCenterFuncClipped_steady;
x1 = ceil(10*max(AeroCenterFuncClipped_AeroCenterFuncs(:)))/10;
y0 = dev_wb_steady_bins_meanCIstd(n_loc);
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('AeroCenterFuncClipped','fontsize',10)
ylabel('Clipped wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([x0 x1 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% coeff
subplot(3,2,4)
plot((AeroCenterFuncClipped_AeroCenterFuncs),devMOD_wb_R_AeroCenterFuncs_bins(n_loc,:),'.b')

RC = devMOD_wb_R_AeroCenterFuncs_bins_meanCIstd(n_loc,1);
plot([x0 x1],[RC RC],'-r','linewidth',2)
xlabel('AeroCenterFuncClipped','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -180;
ang_max = 540;
axis([x0 x1 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% left - right wing = intact - clipped
subplot(3,2,5)
plot((AeroCenterFuncRatio),Ddev_wb_bins(n_loc,:),'.g')
plot((AeroCenterFuncRatio_AeroCenterFuncs),Ddev_wb_AeroCenterFuncs_bins(n_loc,:),'.b')

% plot trend
RC = DdevMOD_wb_AeroCenterFuncs_bins_meanCIstd(n_loc,1);
x0 = AeroCenterFuncRatio_steady;
x1 = ceil(10*max(AeroCenterFuncRatio_AeroCenterFuncs(:)))/10;
y0 = 0;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('AeroCenterFuncRatio','fontsize',10)
ylabel('R - L','fontsize',10)

ang_min = -60;
ang_max = 30;
axis([x0 x1 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% coeff
subplot(3,2,6)
plot((AeroCenterFuncRatio_AeroCenterFuncs),DdevMOD_wb_AeroCenterFuncs_bins(n_loc,:),'.b')

RC = DdevMOD_wb_AeroCenterFuncs_bins_meanCIstd(n_loc,1);
plot([x0 x1],[RC RC],'-r','linewidth',2)
xlabel('AeroCenterFuncRatio','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -270;
ang_max = 90;
axis([x0 x1 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max)
