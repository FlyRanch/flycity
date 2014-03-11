
figure
subplot(2,2,1)
hold on
subplot(2,2,2)
hold on
subplot(2,2,3)
hold on
subplot(2,2,4)
hold on

%% stroke max
subplot(2,2,1)
plot(rot_dot_dot_R_mean_wb,stroke_wb_L_bins(1,:),'.k')
plot(rot_dot_dot_R_mean_wb,stroke_wb_R_bins(1,:),'.k')
plot(rotRaccel_rotAxisRaccel,stroke_wb_L_rotAxisRaccel_bins(1,:),'.b')
plot(rotRaccel_rotAxisRaccel,stroke_wb_R_rotAxisRaccel_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_rotAxisRaccel_bins_meanCIstd(1,1)/rotRaccel_norm;
var0 = stroke_wb_steady_bins_meanCIstd(1,1);
x0 = 0;
x1 = max(rotRaccel_rotAxisRaccel);
var1 = var0 + RC * (x1-x0);
xmin1 = min(rotRaccel_rotAxisRaccel);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Accel','fontsize',10)
ylabel('stroke max','fontsize',10)

ang_min = 45;
ang_max = 135;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,3)
plot(rotRaccel_rotAxisRaccel,strokeMOD_wb_L_rotAxisRaccel_bins(1,:),'.b')
plot(rotRaccel_rotAxisRaccel,strokeMOD_wb_R_rotAxisRaccel_bins(1,:),'.b')

% plot mean
RC = strokeMOD_wb_rotAxisRaccel_bins_meanCIstd(1,1);
x0 = min(rotRaccel_rotAxisRaccel); 
x1 = max(rotRaccel_rotAxisRaccel); 
var0 = RC;
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('rotR Accel','fontsize',10)
ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% stroke min
n_loc = round(Rds_steady_meanCIstd(1)*size(pitch_wb_L_bins,1));

subplot(2,2,2)
plot(rot_dot_dot_R_mean_wb,stroke_wb_L_bins(n_loc,:),'.k')
plot(rot_dot_dot_R_mean_wb,stroke_wb_R_bins(n_loc,:),'.k')
plot(rotRaccel_rotAxisRaccel,stroke_wb_L_rotAxisRaccel_bins(n_loc,:),'.b')
plot(rotRaccel_rotAxisRaccel,stroke_wb_R_rotAxisRaccel_bins(n_loc,:),'.b')

% plot trend
RC = strokeMOD_wb_rotAxisRaccel_bins_meanCIstd(n_loc,1)/rotRaccel_norm;
var0 = stroke_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0; 
x1 = max(rotRaccel_rotAxisRaccel); 
var1 = var0 + RC * (x1-x0);
xmin1 = min(rotRaccel_rotAxisRaccel); 
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Accel','fontsize',10)
ylabel('stroke min','fontsize',10)

ang_min = -90;
ang_max = 0;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,4)
plot(rotRaccel_rotAxisRaccel,strokeMOD_wb_L_rotAxisRaccel_bins(n_loc,:),'.b')
plot(rotRaccel_rotAxisRaccel,strokeMOD_wb_R_rotAxisRaccel_bins(n_loc,:),'.b')

% plot mean
RC = strokeMOD_wb_rotAxisRaccel_bins_meanCIstd(n_loc,1);
x0 = min(rotRaccel_rotAxisRaccel); 
x1 = max(rotRaccel_rotAxisRaccel); 
var0 = RC;
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('rotR Accel','fontsize',10)
ylabel('strokemin mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

