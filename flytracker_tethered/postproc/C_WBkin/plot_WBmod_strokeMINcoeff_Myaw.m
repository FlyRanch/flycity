%% stroke min

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
% plot(abs(Myaw_mean_wb),stroke_wb_L_bins(n_loc,:),'.k')
% plot(abs(Myaw_mean_wb),stroke_wb_R_bins(n_loc,:),'.k')
plot(abs(Myaw_YawAccel),stroke_wb_fwd_YawAccel_bins(n_loc,:),'.b')

% plot trend
RC = strokeMOD_wb_fwd_YawAccel_bins_meanCIstd(n_loc,1)/Myaw_norm;
x0 = 0;
var0 = stroke_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max(abs(Myaw_YawAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Acceleration','fontsize',10)
ylabel('forwards wing','fontsize',10)

ang_min = -90;
ang_max = 0;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
% plot(abs(Myaw_mean_wb),stroke_wb_L_bins(n_loc,:),'.k')
% plot(abs(Myaw_mean_wb),stroke_wb_R_bins(n_loc,:),'.k')
plot(abs(Myaw_YawAccel),stroke_wb_rwd_YawAccel_bins(n_loc,:),'.b')

% plot trend
RC = strokeMOD_wb_rwd_YawAccel_bins_meanCIstd(n_loc,1)/Myaw_norm;
x0 = 0;
var0 = stroke_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max(abs(Myaw_YawAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Acceleration','fontsize',10)
ylabel('rearward wing','fontsize',10)

ang_min = -90;
ang_max = 0;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
% plot(abs(Myaw_mean_wb),stroke_wb_L_bins(n_loc,:),'.k')
% plot(abs(Myaw_mean_wb),stroke_wb_R_bins(n_loc,:),'.k')
plot(abs(Myaw_YawAccel),Dstroke_wb_YawAccel_bins(n_loc,:),'.b')

% plot trend
RC = DstrokeMOD_wb_YawAccel_bins_meanCIstd(n_loc,1)/Myaw_norm;
x0 = 0;
var0 = 0;
x1 = max(abs(Myaw_YawAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Yaw Acceleration','fontsize',10)
ylabel('fwd - rwd','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% coeff
% for
subplot(3,2,2)
plot(abs(Myaw_YawAccel),strokeMOD_wb_fwd_YawAccel_bins(n_loc,:),'.b')

RC = strokeMOD_wb_fwd_YawAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(Myaw_YawAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Acceleration','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% rear
subplot(3,2,4)
plot(abs(Myaw_YawAccel),strokeMOD_wb_rwd_YawAccel_bins(n_loc,:),'.b')

RC = strokeMOD_wb_rwd_YawAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(Myaw_YawAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Acceleration','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% for rear
subplot(3,2,6)
plot(abs(Myaw_YawAccel),DstrokeMOD_wb_YawAccel_bins(n_loc,:),'.b')

RC = DstrokeMOD_wb_YawAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(Myaw_YawAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Yaw Acceleration','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


