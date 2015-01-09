%% dev max

figure
subplot(3,2,1)
title('devmax@upstroke')
hold on
subplot(3,2,2)
title('devmax@up mod coeff')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

start = 50;
stop = 150;
var_loc = dev_wb_steady_bins_meanCIstd(start:stop,1);
n_loc = start+find(var_loc == max(var_loc));

subplot(3,2,1)
% plot((rot_dot_dot_R_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
plot((rot_dot_dot_R_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot(rotRaccel_rotAxisRaccel,dev_wb_R_rotAxisRaccel_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_R_rotAxisRaccel_bins_meanCIstd(n_loc,1)/rotRaccel_norm;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0;
x1 = max(rotRaccel_rotAxisRaccel);
var1 = var0 + RC * (x1-x0);
xmin1 = min(rotRaccel_rotAxisRaccel);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Acceleration','fontsize',10)
ylabel('right wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
plot((rot_dot_dot_R_mean_wb),dev_wb_L_bins(n_loc,:),'.k')
% plot((rot_dot_dot_R_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot((rotRaccel_rotAxisRaccel),dev_wb_L_rotAxisRaccel_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_L_rotAxisRaccel_bins_meanCIstd(n_loc,1)/rotRaccel_norm;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0;
x1 = max(rotRaccel_rotAxisRaccel);
var1 = var0 + RC * (x1-x0);
xmin1 = min(rotRaccel_rotAxisRaccel);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Acceleration','fontsize',10)
ylabel('left wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
plot((rot_dot_dot_R_mean_wb),Ddev_wb_bins(n_loc,:),'.k')
% plot((rot_dot_dot_R_mean_wb),dev_wb_R_bins(n_loc,:),'.k')
plot((rotRaccel_rotAxisRaccel),Ddev_wb_rotAxisRaccel_bins(n_loc,:),'.b')

% plot trend
RC = DdevMOD_wb_rotAxisRaccel_bins_meanCIstd(n_loc,1)/rotRaccel_norm;
var0 = 0;
x0 = 0;
x1 = max(rotRaccel_rotAxisRaccel);
var1 = var0 + RC * (x1-x0);
xmin1 = min(rotRaccel_rotAxisRaccel);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
xlabel('rotR Acceleration','fontsize',10)
ylabel('R - L','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% coeff
% right
subplot(3,2,2)
plot((rotRaccel_rotAxisRaccel),devMOD_wb_R_rotAxisRaccel_bins(n_loc,:),'.b')

RC = devMOD_wb_R_rotAxisRaccel_bins_meanCIstd(n_loc,1);
var0 = 0;
x0 = 0; 
x1 = max(rotRaccel_rotAxisRaccel); 
var1 = RC
xmin1 = min(rotRaccel_rotAxisRaccel); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Acceleration','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Left
subplot(3,2,4)
plot((rotRaccel_rotAxisRaccel),devMOD_wb_L_rotAxisRaccel_bins(n_loc,:),'.b')

RC = devMOD_wb_L_rotAxisRaccel_bins_meanCIstd(n_loc,1);
var0 = 0;
x0 = 0; 
x1 = max(rotRaccel_rotAxisRaccel); 
var1 = RC
xmin1 = min(rotRaccel_rotAxisRaccel); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Acceleration','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% right-Left
subplot(3,2,6)
plot((rotRaccel_rotAxisRaccel),DdevMOD_wb_rotAxisRaccel_bins(n_loc,:),'.b')

RC = DdevMOD_wb_rotAxisRaccel_bins_meanCIstd(n_loc,1);
var0 = 0;
x0 = 0; 
x1 = max(rotRaccel_rotAxisRaccel); 
var1 = RC
xmin1 = min(rotRaccel_rotAxisRaccel); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
xlabel('rotR Acceleration','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


