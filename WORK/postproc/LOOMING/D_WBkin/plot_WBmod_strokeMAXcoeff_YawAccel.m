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

subplot(3,2,1)
% plot(abs(yaw_dot_dot_mean_wb),stroke_wb_L_bins(1,:),'.k')
% plot(abs(yaw_dot_dot_mean_wb),stroke_wb_R_bins(1,:),'.k')
plot(abs(yawaccel_YawAccel),stroke_wb_fwd_YawAccel_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_fwd_YawAccel_bins_meanCIstd(1,1)/yawaccel_norm;
x0 = 0;
var0 = stroke_wb_steady_bins_meanCIstd(1,1);
x1 = max(abs(yawaccel_YawAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Acceleration','fontsize',10)
ylabel('forwards wing','fontsize',10)

ang_min = 45;
ang_max = 135;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
% plot(abs(yaw_dot_dot_mean_wb),stroke_wb_L_bins(1,:),'.k')
% plot(abs(yaw_dot_dot_mean_wb),stroke_wb_R_bins(1,:),'.k')
plot(abs(yawaccel_YawAccel),stroke_wb_rwd_YawAccel_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_rwd_YawAccel_bins_meanCIstd(1,1)/yawaccel_norm;
x0 = 0;
var0 = stroke_wb_steady_bins_meanCIstd(1,1);
x1 = max(abs(yawaccel_YawAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Acceleration','fontsize',10)
ylabel('rearward wing','fontsize',10)

ang_min = 45;
ang_max = 135;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
% plot(abs(yaw_dot_dot_mean_wb),stroke_wb_L_bins(1,:),'.k')
% plot(abs(yaw_dot_dot_mean_wb),stroke_wb_R_bins(1,:),'.k')
plot(abs(yawaccel_YawAccel),Dstroke_wb_YawAccel_bins(1,:),'.b')

% plot trend
RC = DstrokeMOD_wb_YawAccel_bins_meanCIstd(1,1)/yawaccel_norm;
x0 = 0;
var0 = 0;
x1 = max(abs(yawaccel_YawAccel));
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
plot(abs(yawaccel_YawAccel),strokeMOD_wb_fwd_YawAccel_bins(1,:),'.b')

RC = strokeMOD_wb_fwd_YawAccel_bins_meanCIstd(1,1);
x0 = 0;
var0 = RC;
x1 = max(abs(yawaccel_YawAccel));
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
plot(abs(yawaccel_YawAccel),strokeMOD_wb_rwd_YawAccel_bins(1,:),'.b')

RC = strokeMOD_wb_rwd_YawAccel_bins_meanCIstd(1,1);
x0 = 0;
var0 = RC;
x1 = max(abs(yawaccel_YawAccel));
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
plot(abs(yawaccel_YawAccel),DstrokeMOD_wb_YawAccel_bins(1,:),'.b')

RC = DstrokeMOD_wb_YawAccel_bins_meanCIstd(1,1);
x0 = 0;
var0 = RC;
x1 = max(abs(yawaccel_YawAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Yaw Acceleration','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


