%% dev max

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

subplot(3,2,1)
% plot((M_axis1_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
plot((M_axis1_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot(M_axis1_TorqueAxis1UP,dev_wb_R_TorqueAxis1UP_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_R_TorqueAxis1UP_bins_meanCIstd(n_loc,1)/M_axis1_norm;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0;
x1 = max(M_axis1_TorqueAxis1UP);
var1 = var0 + RC * (x1-x0);
xmin1 = min(M_axis1_TorqueAxis1UP);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
ylabel('right wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
plot((M_axis1_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
% plot((M_axis1_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot((M_axis1_TorqueAxis1UP),dev_wb_L_TorqueAxis1UP_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_L_TorqueAxis1UP_bins_meanCIstd(n_loc,1)/M_axis1_norm;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0;
x1 = max(M_axis1_TorqueAxis1UP);
var1 = var0 + RC * (x1-x0);
xmin1 = min(M_axis1_TorqueAxis1UP);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
ylabel('left wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
plot((M_axis1_mean_wb),Ddev_wb_bins(n_loc,:),'.k')
% plot((M_axis1_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot((M_axis1_TorqueAxis1UP),Ddev_wb_TorqueAxis1UP_bins(n_loc,:),'.b')

% plot trend
RC = DdevMOD_wb_TorqueAxis1UP_bins_meanCIstd(n_loc,1)/M_axis1_norm;
var0 = 0;
x0 = 0;
x1 = max(M_axis1_TorqueAxis1UP);
var1 = var0 + RC * (x1-x0);
xmin1 = min(M_axis1_TorqueAxis1UP);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
xlabel('rotR Torque','fontsize',10)
ylabel('R - L','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% coeff
% right
subplot(3,2,2)
plot((M_axis1_TorqueAxis1UP),devMOD_wb_R_TorqueAxis1UP_bins(n_loc,:),'.b')

RC = devMOD_wb_R_TorqueAxis1UP_bins_meanCIstd(n_loc,1);
var0 = 0;
x0 = 0; 
x1 = max(M_axis1_TorqueAxis1UP); 
var1 = RC
xmin1 = min(M_axis1_TorqueAxis1UP); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Left
subplot(3,2,4)
plot((M_axis1_TorqueAxis1UP),devMOD_wb_L_TorqueAxis1UP_bins(n_loc,:),'.b')

RC = devMOD_wb_L_TorqueAxis1UP_bins_meanCIstd(n_loc,1);
var0 = 0;
x0 = 0; 
x1 = max(M_axis1_TorqueAxis1UP); 
var1 = RC
xmin1 = min(M_axis1_TorqueAxis1UP); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% right-Left
subplot(3,2,6)
plot((M_axis1_TorqueAxis1UP),DdevMOD_wb_TorqueAxis1UP_bins(n_loc,:),'.b')

RC = DdevMOD_wb_TorqueAxis1UP_bins_meanCIstd(n_loc,1);
var0 = 0;
x0 = 0; 
x1 = max(M_axis1_TorqueAxis1UP); 
var1 = RC
xmin1 = min(M_axis1_TorqueAxis1UP); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
xlabel('rotR Torque','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


