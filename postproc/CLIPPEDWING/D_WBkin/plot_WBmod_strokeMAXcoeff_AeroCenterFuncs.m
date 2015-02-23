%% stroke max

figure
subplot(3,2,1)
title('stroke max')
hold on
subplot(3,2,2)
title('strokemax mod coefficient')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

%% left wing = intact
subplot(3,2,1)
plot((AeroCenterFuncIntact),stroke_wb_L_bins(1,:),'.g')
plot((AeroCenterFuncIntact_AeroCenterFuncs),stroke_wb_L_AeroCenterFuncs_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_L_AeroCenterFuncs_bins_meanCIstd(1,1);
x0 = AeroCenterFuncIntact_steady;
x1 = floor(10*min(AeroCenterFuncIntact_AeroCenterFuncs(:)))/10;
y0 = stroke_wb_steady_bins_meanCIstd(1);
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('AeroCenterFuncIntact','fontsize',10)
ylabel('Intact wing','fontsize',10)

ang_min = 45;
ang_max = 90;
axis([x1 x0 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% coeff
subplot(3,2,2)
plot((AeroCenterFuncIntact_AeroCenterFuncs),strokeMOD_wb_L_AeroCenterFuncs_bins(1,:),'.b')

RC = strokeMOD_wb_L_AeroCenterFuncs_bins_meanCIstd(1,1);
plot([x0 x1],[RC RC],'-r','linewidth',2)
xlabel('AeroCenterFuncIntact','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -90;
ang_max = 90;
axis([x1 x0 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% right wing = clipped
subplot(3,2,3)
plot((AeroCenterFuncClipped),stroke_wb_R_bins(1,:),'.g')
plot(AeroCenterFuncClipped_AeroCenterFuncs,stroke_wb_R_AeroCenterFuncs_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_R_AeroCenterFuncs_bins_meanCIstd(1,1);
x0 = AeroCenterFuncClipped_steady;
x1 = ceil(10*max(AeroCenterFuncClipped_AeroCenterFuncs(:)))/10;
y0 = stroke_wb_steady_bins_meanCIstd(1);
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('AeroCenterFuncClipped','fontsize',10)
ylabel('Clipped wing','fontsize',10)

ang_min = 45;
ang_max = 90;
axis([x0 x1 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% coeff
subplot(3,2,4)
plot((AeroCenterFuncClipped_AeroCenterFuncs),strokeMOD_wb_R_AeroCenterFuncs_bins(1,:),'.b')

RC = strokeMOD_wb_R_AeroCenterFuncs_bins_meanCIstd(1,1);
plot([x0 x1],[RC RC],'-r','linewidth',2)
xlabel('AeroCenterFuncClipped','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -90;
ang_max = 90;
axis([x0 x1 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% left - right wing = intact - clipped
subplot(3,2,5)
plot((AeroCenterFuncRatio),Dstroke_wb_bins(1,:),'.g')
plot((AeroCenterFuncRatio_AeroCenterFuncs),Dstroke_wb_AeroCenterFuncs_bins(1,:),'.b')

% plot trend
RC = DstrokeMOD_wb_AeroCenterFuncs_bins_meanCIstd(1,1);
x0 = AeroCenterFuncRatio_steady;
x1 = ceil(10*max(AeroCenterFuncRatio_AeroCenterFuncs(:)))/10;
y0 = 0;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('AeroCenterFuncRatio','fontsize',10)
ylabel('R - L','fontsize',10)

ang_min = -15;
ang_max = 15;
axis([x0 x1 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% coeff
subplot(3,2,6)
plot((AeroCenterFuncRatio_AeroCenterFuncs),DstrokeMOD_wb_AeroCenterFuncs_bins(1,:),'.b')

RC = DstrokeMOD_wb_AeroCenterFuncs_bins_meanCIstd(1,1);
plot([x0 x1],[RC RC],'-r','linewidth',2)
xlabel('AeroCenterFuncRatio','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -30;
ang_max = 30;
axis([x0 x1 ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max)
