%% dev min up

figure
subplot(3,2,1)
title('devmin@upstroke')
hold on
subplot(3,2,2)
title('devmin@upstroke mod coeff')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

% min up dev bin location
start = 100;
stop = 200;
var_loc = dev_wb_steady_bins_meanCIstd(start:stop,1);
n_loc = start + find(var_loc == min(var_loc));

% upwards
subplot(3,2,1)
% plot(abs(roll_dot_dot_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
% plot(abs(roll_dot_dot_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot(abs(rollaccel_RollAccel),dev_wb_up_RollAccel_bins(n_loc,:),'.b')

RC = devMOD_wb_up_RollAccel_bins_meanCIstd(n_loc,1)/rollaccel_norm;
x0 = 0;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max(abs(rollaccel_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
ylabel('upwards wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% downwards
subplot(3,2,3)
% plot(abs(roll_dot_dot_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
% plot(abs(roll_dot_dot_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot(abs(rollaccel_RollAccel),dev_wb_down_RollAccel_bins(n_loc,:),'.b')

RC = devMOD_wb_down_RollAccel_bins_meanCIstd(n_loc,1)/rollaccel_norm;
x0 = 0;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max(abs(rollaccel_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
ylabel('downward wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% up - down
subplot(3,2,5)
% plot(abs(roll_dot_dot_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
% plot(abs(roll_dot_dot_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot(abs(rollaccel_RollAccel),Ddev_wb_RollAccel_bins(n_loc,:),'.b')

RC = DdevMOD_wb_RollAccel_bins_meanCIstd(n_loc,1)/rollaccel_norm;
x0 = 0;
var0 = 0;
x1 = max(abs(rollaccel_RollAccel));
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
plot(abs(rollaccel_RollAccel),devMOD_wb_up_RollAccel_bins(n_loc,:),'.b')

RC = devMOD_wb_up_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(rollaccel_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
% ylabel('devmin mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% down
subplot(3,2,4)
plot(abs(rollaccel_RollAccel),devMOD_wb_down_RollAccel_bins(n_loc,:),'.b')

RC = devMOD_wb_down_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(rollaccel_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
% ylabel('devmin mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% up - down
subplot(3,2,6)
plot(abs(rollaccel_RollAccel),DdevMOD_wb_RollAccel_bins(n_loc,:),'.b')

RC = DdevMOD_wb_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(rollaccel_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Roll Acceleration','fontsize',10)
% ylabel('devmin mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

