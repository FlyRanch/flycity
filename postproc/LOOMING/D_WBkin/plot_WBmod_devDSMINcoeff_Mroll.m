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
% plot(abs(Mroll_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
% plot(abs(Mroll_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot(abs(Mroll_RollAccel),dev_wb_up_RollAccel_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_up_RollAccel_bins_meanCIstd(n_loc,1)/Mroll_norm;
x0 = 0;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max(abs(Mroll_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
ylabel('upwards wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
% plot(abs(Mroll_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
% plot(abs(Mroll_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot(abs(Mroll_RollAccel),dev_wb_down_RollAccel_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_down_RollAccel_bins_meanCIstd(n_loc,1)/Mroll_norm;
x0 = 0;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = max(abs(Mroll_RollAccel));
var1 = var0 + RC * (x1-x0);
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
ylabel('downward wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
% plot(abs(Mroll_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
% plot(abs(Mroll_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot(abs(Mroll_RollAccel),Ddev_wb_RollAccel_bins(n_loc,:),'.b')

% plot trend
RC = DdevMOD_wb_RollAccel_bins_meanCIstd(n_loc,1)/Mroll_norm;
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
plot(abs(Mroll_RollAccel),devMOD_wb_up_RollAccel_bins(n_loc,:),'.b')

RC = devMOD_wb_up_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(Mroll_RollAccel));
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
plot(abs(Mroll_RollAccel),devMOD_wb_down_RollAccel_bins(n_loc,:),'.b')

RC = devMOD_wb_down_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(Mroll_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
% xlabel('Roll Acceleration','fontsize',10)
% ylabel('devmin mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% up down
subplot(3,2,6)
plot(abs(Mroll_RollAccel),DdevMOD_wb_RollAccel_bins(n_loc,:),'.b')

RC = DdevMOD_wb_RollAccel_bins_meanCIstd(n_loc,1);
x0 = 0;
var0 = RC;
x1 = max(abs(Mroll_RollAccel));
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Roll Acceleration','fontsize',10)
% ylabel('devmin mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


