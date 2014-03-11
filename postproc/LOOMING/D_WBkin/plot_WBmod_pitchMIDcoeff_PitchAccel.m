
figure
subplot(2,2,1)
hold on
subplot(2,2,2)
hold on
subplot(2,2,3)
hold on
subplot(2,2,4)
hold on

%% pitch mid downstroke
subplot(2,2,1)
n_loc = round(Rds_steady_meanCIstd(1) * size(pitch_wb_L_bins,1)/2);
plot(pitch_dot_dot_mean_wb,pitch_wb_L_bins(n_loc,:),'.k')
plot(pitch_dot_dot_mean_wb,pitch_wb_R_bins(n_loc,:),'.k')
plot(pitchaccel_PitchAccel,pitch_wb_L_PitchAccel_bins(n_loc,:),'.b')
plot(pitchaccel_PitchAccel,pitch_wb_R_PitchAccel_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_PitchAccel_bins_meanCIstd(n_loc,1)/pitchaccel_norm;
var0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0;
x1 = max(pitchaccel_PitchAccel);
var1 = var0 + RC * (x1-x0);
xmin1 = min(pitchaccel_PitchAccel);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('Pitch Accel','fontsize',10)
ylabel('pitch mid down','fontsize',10)

ang_min = 90;
ang_max = 180;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,3)
plot(pitchaccel_PitchAccel,pitchMOD_wb_L_PitchAccel_bins(n_loc,:),'.b')
plot(pitchaccel_PitchAccel,pitchMOD_wb_R_PitchAccel_bins(n_loc,:),'.b')

% plot mean
RC = pitchMOD_wb_PitchAccel_bins_meanCIstd(n_loc,1);
x0 = min(pitchaccel_PitchAccel); 
x1 = max(pitchaccel_PitchAccel); 
var0 = RC;
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Pitch Accel','fontsize',10)
ylabel('pitch mid down-F mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% pitch mid upstroke
subplot(2,2,2)
n_loc = round((Rds_steady_meanCIstd(1)+1)*size(pitch_wb_L_bins,1)/2);
plot(pitch_dot_dot_mean_wb,pitch_wb_L_bins(n_loc,:),'.k')
plot(pitch_dot_dot_mean_wb,pitch_wb_R_bins(n_loc,:),'.k')
plot(pitchaccel_PitchAccel,pitch_wb_L_PitchAccel_bins(n_loc,:),'.b')
plot(pitchaccel_PitchAccel,pitch_wb_R_PitchAccel_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_PitchAccel_bins_meanCIstd(n_loc,1)/pitchaccel_norm;
var0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0; 
x1 = max(pitchaccel_PitchAccel); 
var1 = var0 + RC * (x1-x0);
xmin1 = min(pitchaccel_PitchAccel); 
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('Pitch Accel','fontsize',10)
ylabel('pitch mid up','fontsize',10)

ang_min = 0;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,4)
plot(pitchaccel_PitchAccel,pitchMOD_wb_L_PitchAccel_bins(n_loc,:),'.b')
plot(pitchaccel_PitchAccel,pitchMOD_wb_R_PitchAccel_bins(n_loc,:),'.b')

% plot mean
RC = pitchMOD_wb_PitchAccel_bins_meanCIstd(n_loc,1);
x0 = min(pitchaccel_PitchAccel); 
x1 = max(pitchaccel_PitchAccel); 
var0 = RC;
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Pitch Accel','fontsize',10)
ylabel('pitch mid up-F mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 