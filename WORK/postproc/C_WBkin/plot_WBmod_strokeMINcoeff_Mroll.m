%% stroke max

figure
subplot(3,2,1)
title('min stroke')
hold on
subplot(3,2,2)
title('stroke min mod coefficient')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

n_loc = round(Rds_steady_meanCIstd(1)*size(pitch_wb_L_bins,1));

subplot(3,2,1)
% plot(abs(Mroll_mean_wb),stroke_wb_L_bins(n_loc,:),'.k')
% plot(abs(Mroll_mean_wb),stroke_wb_R_bins(n_loc,:),'.k')
plot(abs(Mroll_RollAccel),stroke_wb_up_RollAccel_bins(n_loc,:),'.b')

% plot trend
RC = strokeMOD_wb_up_RollAccel_bins_meanCIstd(n_loc,1)/Mroll_norm;
x0 = 0;
var0 = stroke_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max(abs(Mroll_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
ylabel('upwards wing','fontsize',10)

ang_min = -90;
ang_max = 0;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
% plot(abs(Mroll_mean_wb),stroke_wb_L_bins(n_loc,:),'.k')
% plot(abs(Mroll_mean_wb),stroke_wb_R_bins(n_loc,:),'.k')
plot(abs(Mroll_RollAccel),stroke_wb_down_RollAccel_bins(n_loc,:),'.b')

% plot trend
RC = strokeMOD_wb_down_RollAccel_bins_meanCIstd(n_loc,1)/Mroll_norm;
x0 = 0;
var0 = stroke_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max(abs(Mroll_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
ylabel('downward wing','fontsize',10)

ang_min = -90;
ang_max = 0;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
% plot(abs(Mroll_mean_wb),stroke_wb_L_bins(n_loc,:),'.k')
% plot(abs(Mroll_mean_wb),stroke_wb_R_bins(n_loc,:),'.k')
plot(abs(Mroll_RollAccel),Dstroke_wb_RollAccel_bins(n_loc,:),'.b')

% plot trend
RC = DstrokeMOD_wb_RollAccel_bins_meanCIstd(n_loc,1)/Mroll_norm;
x0 = 0;
var0 = 0;
x1 = max(abs(Mroll_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Roll Acceleration','fontsize',10)
ylabel('up - down','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% coeff
% up
subplot(3,2,2)
plot(abs(Mroll_RollAccel),strokeMOD_wb_up_RollAccel_bins(n_loc,:),'.b')

RC = strokeMOD_wb_up_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(Mroll_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% down
subplot(3,2,4)
plot(abs(Mroll_RollAccel),strokeMOD_wb_down_RollAccel_bins(n_loc,:),'.b')

RC = strokeMOD_wb_down_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(Mroll_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% up down
subplot(3,2,6)
plot(abs(Mroll_RollAccel),DstrokeMOD_wb_RollAccel_bins(n_loc,:),'.b')

RC = DstrokeMOD_wb_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(Mroll_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Roll Acceleration','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


