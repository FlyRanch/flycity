
figure
subplot(2,2,1)
hold on
subplot(2,2,2)
hold on
subplot(2,2,3)
hold on
subplot(2,2,4)
hold on

%% max dev upstroke
start = 50;
stop = 150;
var_loc = dev_wb_steady_bins_meanCIstd(start:stop,1);
n_loc = start+find(var_loc == max(var_loc));

subplot(2,2,1)
plot(pitch_dot_dot_mean_wb,dev_wb_L_bins(n_loc,:),'.k')
plot(pitch_dot_dot_mean_wb,dev_wb_R_bins(n_loc,:),'.k')
plot(pitchaccel_PitchAccel,dev_wb_L_PitchAccel_bins(n_loc,:),'.b')
plot(pitchaccel_PitchAccel,dev_wb_R_PitchAccel_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_PitchAccel_bins_meanCIstd(n_loc,1)/pitchaccel_norm;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0;
x1 = max(pitchaccel_PitchAccel);
var1 = var0 + RC * (x1-x0);
xmin1 = min(pitchaccel_PitchAccel);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('Pitch Accel','fontsize',10)
ylabel('max dev up','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,3)
plot(pitchaccel_PitchAccel,devMOD_wb_L_PitchAccel_bins(n_loc,:),'.b')
plot(pitchaccel_PitchAccel,devMOD_wb_R_PitchAccel_bins(n_loc,:),'.b')

% plot mean
RC = devMOD_wb_PitchAccel_bins_meanCIstd(n_loc,1);
x0 = min(pitchaccel_PitchAccel); 
x1 = max(pitchaccel_PitchAccel); 
var0 = RC;
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Pitch Accel','fontsize',10)
ylabel('max dev up mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% min dev up
start = 100;
stop = 200;
var_loc = dev_wb_steady_bins_meanCIstd(start:stop,1);
n_loc = start + find(var_loc == min(var_loc));

subplot(2,2,2)
plot(pitch_dot_dot_mean_wb,dev_wb_L_bins(n_loc,:),'.k')
plot(pitch_dot_dot_mean_wb,dev_wb_R_bins(n_loc,:),'.k')
plot(pitchaccel_PitchAccel,dev_wb_L_PitchAccel_bins(n_loc,:),'.b')
plot(pitchaccel_PitchAccel,dev_wb_R_PitchAccel_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_PitchAccel_bins_meanCIstd(n_loc,1)/pitchaccel_norm;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0;
x1 = max(pitchaccel_PitchAccel);
var1 = var0 + RC * (x1-x0);
xmin1 = min(pitchaccel_PitchAccel);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('Pitch Accel','fontsize',10)
ylabel('min dev up','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,4)
plot(pitchaccel_PitchAccel,devMOD_wb_L_PitchAccel_bins(n_loc,:),'.b')
plot(pitchaccel_PitchAccel,devMOD_wb_R_PitchAccel_bins(n_loc,:),'.b')

% plot mean
RC = devMOD_wb_PitchAccel_bins_meanCIstd(n_loc,1);
x0 = min(pitchaccel_PitchAccel); 
x1 = max(pitchaccel_PitchAccel); 
var0 = RC;
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Pitch Accel','fontsize',10)
ylabel('min dev up-F mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

