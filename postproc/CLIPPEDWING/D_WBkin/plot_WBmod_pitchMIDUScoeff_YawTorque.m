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
n_loc = round((Rds_steady_meanCIstd(1)+1)*size(pitch_wb_fwd_YawTorque_bins,1)/2);

subplot(3,2,1)
plot((Myaw_mean_wb),pitch_wb_L_bins(n_loc,:),'.k')
% plot((Myaw_mean_wb),pitch_wb_R_bins(n_loc,:),'.k')
plot((Myaw_YawTorque),pitch_wb_fwd_YawTorque_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_fwd_YawTorque_bins_meanCIstd(n_loc,1)/Myaw_norm;
x0 = 0;
var0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max((Myaw_YawTorque));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Torque','fontsize',10)
ylabel('forwards wing','fontsize',10)

ang_min = 0;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
% plot((Myaw_mean_wb),pitch_wb_L_bins(n_loc,:),'.k')
plot((Myaw_mean_wb),pitch_wb_R_bins(n_loc,:),'.k')
plot((Myaw_YawTorque),pitch_wb_rwd_YawTorque_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_rwd_YawTorque_bins_meanCIstd(n_loc,1)/Myaw_norm;
x0 = 0;
var0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max((Myaw_YawTorque));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Torque','fontsize',10)
ylabel('rearward wing','fontsize',10)

ang_min = 0;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
% plot((Myaw_mean_wb),pitch_wb_L_bins(n_loc,:),'.k')
plot((Myaw_mean_wb),Dpitch_wb_bins(n_loc,:),'.k')
plot((Myaw_YawTorque),Dpitch_wb_YawTorque_bins(n_loc,:),'.b')

% plot trend
RC = DpitchMOD_wb_YawTorque_bins_meanCIstd(n_loc,1)/Myaw_norm;
x0 = 0;
var0 = 0;
x1 = max((Myaw_YawTorque));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Yaw Torque','fontsize',10)
ylabel('fwd - rwd','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% coeff
% fwd
subplot(3,2,2)
plot((Myaw_YawTorque),pitchMOD_wb_fwd_YawTorque_bins(n_loc,:),'.b')

RC = pitchMOD_wb_fwd_YawTorque_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max((Myaw_YawTorque));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Torque','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% rwd
subplot(3,2,4)
plot((Myaw_YawTorque),pitchMOD_wb_rwd_YawTorque_bins(n_loc,:),'.b')

RC = pitchMOD_wb_rwd_YawTorque_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max((Myaw_YawTorque));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Torque','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% fwd - rwd
subplot(3,2,6)
plot((Myaw_YawTorque),DpitchMOD_wb_YawTorque_bins(n_loc,:),'.b')

RC = DpitchMOD_wb_YawTorque_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max((Myaw_YawTorque));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Yaw Torque','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


