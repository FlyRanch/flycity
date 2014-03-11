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
% plot(abs(roll_dot_dot_mean_wb),stroke_wb_L_bins(1,:),'.k')
% plot(abs(roll_dot_dot_mean_wb),stroke_wb_R_bins(1,:),'.k')
plot(abs(rollaccel_RollAccel),stroke_wb_up_RollAccel_bins(1,:),'ob')

% plot trend
RC = strokeMOD_wb_up_RollAccel_bins_meanCIstd(1,1)/rollaccel_norm;
x0 = 0;
var0 = stroke_wb_steady_bins_meanCIstd(1,1);
x1 = max(abs(rollaccel_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Roll Acceleration','fontsize',10)
ylabel('upwards wing','fontsize',10)

subplot(3,2,3)
% plot(abs(roll_dot_dot_mean_wb),stroke_wb_L_bins(1,:),'.k')
% plot(abs(roll_dot_dot_mean_wb),stroke_wb_R_bins(1,:),'.k')
plot(abs(rollaccel_RollAccel),stroke_wb_down_RollAccel_bins(1,:),'ob')

% plot trend
RC = strokeMOD_wb_down_RollAccel_bins_meanCIstd(1,1)/rollaccel_norm;
x0 = 0;
var0 = stroke_wb_steady_bins_meanCIstd(1,1);
x1 = max(abs(rollaccel_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Roll Acceleration','fontsize',10)
ylabel('downward wing','fontsize',10)

subplot(3,2,5)
% plot(abs(roll_dot_dot_mean_wb),stroke_wb_L_bins(1,:),'.k')
% plot(abs(roll_dot_dot_mean_wb),stroke_wb_R_bins(1,:),'.k')
plot(abs(rollaccel_RollAccel),Dstroke_wb_RollAccel_bins(1,:),'ob')

% plot trend
RC = DstrokeMOD_wb_RollAccel_bins_meanCIstd(1,1)/rollaccel_norm;
x0 = 0;
var0 = 0;
x1 = max(abs(rollaccel_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Roll Acceleration','fontsize',10)
ylabel('up - down','fontsize',10)

%% coeff
% up
subplot(3,2,2)
plot(abs(rollaccel_RollAccel),strokeMOD_wb_up_RollAccel_bins(1,:),'ob')

RC = strokeMOD_wb_up_RollAccel_bins_meanCIstd(1,1);
x0 = 0;
var0 = RC;
x1 = max(abs(rollaccel_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Roll Acceleration','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

% down
subplot(3,2,4)
plot(abs(rollaccel_RollAccel),strokeMOD_wb_down_RollAccel_bins(1,:),'ob')

RC = strokeMOD_wb_down_RollAccel_bins_meanCIstd(1,1);
x0 = 0;
var0 = RC;
x1 = max(abs(rollaccel_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Roll Acceleration','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

% up down
subplot(3,2,6)
plot(abs(rollaccel_RollAccel),DstrokeMOD_wb_RollAccel_bins(1,:),'ob')

RC = DstrokeMOD_wb_RollAccel_bins_meanCIstd(1,1);
x0 = 0;
var0 = RC;
x1 = max(abs(rollaccel_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Roll Acceleration','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)


