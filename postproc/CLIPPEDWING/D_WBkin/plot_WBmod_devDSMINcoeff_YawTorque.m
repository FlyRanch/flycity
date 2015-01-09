%% dev min down

figure
subplot(3,2,1)
title('devmin@downstroke')
hold on
subplot(3,2,2)
title('devmin@down mod coeff')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

% min down dev bin location
var_loc = dev_wb_steady_bins_meanCIstd(1:100,1);
n_loc = find(var_loc == min(var_loc));

subplot(3,2,1)
plot((Myaw_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
% plot((Myaw_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot((Myaw_YawTorque),dev_wb_fwd_YawTorque_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_fwd_YawTorque_bins_meanCIstd(n_loc,1)/Myaw_norm;
x0 = 0;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max((Myaw_YawTorque));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Torque','fontsize',10)
ylabel('forwards wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
% plot((Myaw_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
plot((Myaw_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot((Myaw_YawTorque),dev_wb_rwd_YawTorque_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_rwd_YawTorque_bins_meanCIstd(n_loc,1)/Myaw_norm;
x0 = 0;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max((Myaw_YawTorque));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Torque','fontsize',10)
ylabel('rearward wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
% plot((Myaw_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
plot((Myaw_mean_wb),Ddev_wb_bins(n_loc,:),'.k')
plot((Myaw_YawTorque),Ddev_wb_YawTorque_bins(n_loc,:),'.b')

% plot trend
RC = DdevMOD_wb_YawTorque_bins_meanCIstd(n_loc,1)/Myaw_norm;
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
plot((Myaw_YawTorque),devMOD_wb_fwd_YawTorque_bins(n_loc,:),'.b')

RC = devMOD_wb_fwd_YawTorque_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max((Myaw_YawTorque));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Torque','fontsize',10)
% ylabel('devmin mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% rwd
subplot(3,2,4)
plot((Myaw_YawTorque),devMOD_wb_rwd_YawTorque_bins(n_loc,:),'.b')

RC = devMOD_wb_rwd_YawTorque_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max((Myaw_YawTorque));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Yaw Torque','fontsize',10)
% ylabel('devmin mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% fwd - rwd
subplot(3,2,6)
plot((Myaw_YawTorque),DdevMOD_wb_YawTorque_bins(n_loc,:),'.b')

RC = DdevMOD_wb_YawTorque_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max((Myaw_YawTorque));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Yaw Torque','fontsize',10)
% ylabel('devmin mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


