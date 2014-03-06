
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
plot(Mpitch_mean_wb,stroke_wb_L_bins(1,:),'.k')
plot(Mpitch_mean_wb,stroke_wb_R_bins(1,:),'.k')
plot(Mpitch_PitchAccel,stroke_wb_L_PitchAccel_bins(1,:),'.b')
plot(Mpitch_PitchAccel,stroke_wb_R_PitchAccel_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_PitchAccel_bins_meanCIstd(1,1)/Mpitch_norm;
var0 = stroke_wb_steady_bins_meanCIstd(1,1);
x0 = 0;
x1 = max(Mpitch_PitchAccel);
var1 = var0 + RC * (x1-x0);
xmin1 = min(Mpitch_PitchAccel);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('Pitch Accel','fontsize',10)
ylabel('stroke max','fontsize',10)

ang_min = 45;
ang_max = 135;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,3)
plot(Mpitch_PitchAccel,strokeMOD_wb_L_PitchAccel_bins(1,:),'.b')
plot(Mpitch_PitchAccel,strokeMOD_wb_R_PitchAccel_bins(1,:),'.b')

% plot mean
RC = strokeMOD_wb_PitchAccel_bins_meanCIstd(1,1);
x0 = min(Mpitch_PitchAccel); 
x1 = max(Mpitch_PitchAccel); 
var0 = RC;
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Pitch Accel','fontsize',10)
ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% stroke min
n_loc = round(Rds_steady_meanCIstd(1)*size(pitch_wb_L_bins,1));

subplot(2,2,2)
plot(Mpitch_mean_wb,stroke_wb_L_bins(n_loc,:),'.k')
plot(Mpitch_mean_wb,stroke_wb_R_bins(n_loc,:),'.k')
plot(Mpitch_PitchAccel,stroke_wb_L_PitchAccel_bins(n_loc,:),'.b')
plot(Mpitch_PitchAccel,stroke_wb_R_PitchAccel_bins(n_loc,:),'.b')

% plot trend
RC = strokeMOD_wb_PitchAccel_bins_meanCIstd(n_loc,1)/Mpitch_norm;
var0 = stroke_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0; 
x1 = max(Mpitch_PitchAccel); 
var1 = var0 + RC * (x1-x0);
xmin1 = min(Mpitch_PitchAccel); 
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('Pitch Accel','fontsize',10)
ylabel('stroke min','fontsize',10)

ang_min = -90;
ang_max = 0;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,4)
plot(Mpitch_PitchAccel,strokeMOD_wb_L_PitchAccel_bins(n_loc,:),'.b')
plot(Mpitch_PitchAccel,strokeMOD_wb_R_PitchAccel_bins(n_loc,:),'.b')

% plot mean
RC = strokeMOD_wb_PitchAccel_bins_meanCIstd(n_loc,1);
x0 = min(Mpitch_PitchAccel); 
x1 = max(Mpitch_PitchAccel); 
var0 = RC;
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Pitch Accel','fontsize',10)
ylabel('strokemin mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

