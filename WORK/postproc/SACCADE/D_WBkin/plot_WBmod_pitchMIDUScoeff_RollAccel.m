%% pitch @ mid upstroke

figure
subplot(3,2,1)
title('pitch@mid upstroke')
hold on
subplot(3,2,2)
title('pitch@midup mod coeff')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

% mid up bin location
n_loc = round((Rds_steady_meanCIstd(1)+1)*size(pitch_wb_up_RollAccel_bins,1)/2);

subplot(3,2,1)
% plot(abs(roll_dot_dot_mean_wb),pitch_wb_L_bins(n_loc,:),'.k')
% plot(abs(roll_dot_dot_mean_wb),pitch_wb_R_bins(n_loc,:),'.k')
plot(abs(rollaccel_RollAccel),pitch_wb_up_RollAccel_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_up_RollAccel_bins_meanCIstd(n_loc,1)/rollaccel_norm;
x0 = 0;
var0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max(abs(rollaccel_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
ylabel('upwards wing','fontsize',10)

ang_min = 0;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
% plot(abs(roll_dot_dot_mean_wb),pitch_wb_L_bins(n_loc,:),'.k')
% plot(abs(roll_dot_dot_mean_wb),pitch_wb_R_bins(n_loc,:),'.k')
plot(abs(rollaccel_RollAccel),pitch_wb_down_RollAccel_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_down_RollAccel_bins_meanCIstd(n_loc,1)/rollaccel_norm;
x0 = 0;
var0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max(abs(rollaccel_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
ylabel('downward wing','fontsize',10)

ang_min = 0;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
% plot(abs(roll_dot_dot_mean_wb),pitch_wb_L_bins(n_loc,:),'.k')
% plot(abs(roll_dot_dot_mean_wb),pitch_wb_R_bins(n_loc,:),'.k')
plot(abs(rollaccel_RollAccel),Dpitch_wb_RollAccel_bins(n_loc,:),'.b')

% plot trend
RC = DpitchMOD_wb_RollAccel_bins_meanCIstd(n_loc,1)/rollaccel_norm;
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
plot(abs(rollaccel_RollAccel),pitchMOD_wb_up_RollAccel_bins(n_loc,:),'.b')

RC = pitchMOD_wb_up_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(rollaccel_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% down
subplot(3,2,4)
plot(abs(rollaccel_RollAccel),pitchMOD_wb_down_RollAccel_bins(n_loc,:),'.b')

RC = pitchMOD_wb_down_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(rollaccel_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% up down
subplot(3,2,6)
plot(abs(rollaccel_RollAccel),DpitchMOD_wb_RollAccel_bins(n_loc,:),'.b')

RC = DpitchMOD_wb_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(rollaccel_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Roll Acceleration','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


