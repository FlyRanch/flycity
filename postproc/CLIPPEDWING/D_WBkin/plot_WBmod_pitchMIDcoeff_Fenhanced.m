
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
plot(F_mean_wb,pitch_wb_L_bins(n_loc,:),'.k')
plot(F_mean_wb,pitch_wb_R_bins(n_loc,:),'.k')
plot(F_Fenhance,pitch_wb_L_Fenhance_bins(n_loc,:),'.b')
plot(F_Fenhance,pitch_wb_R_Fenhance_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_Fenhance_bins_meanCIstd(n_loc,1)/Fenhance_norm;
F0 = 1;
var0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
F1 = 2;
var1 = var0 + RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
% xlabel('Flight Force','fontsize',10)
ylabel('pitch mid down','fontsize',10)

ang_min = 90;
ang_max = 180;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,3)
plot(F_Fenhance,pitchMOD_wb_L_Fenhance_bins(n_loc,:),'.b')
plot(F_Fenhance,pitchMOD_wb_R_Fenhance_bins(n_loc,:),'.b')

% plot mean
RC = pitchMOD_wb_Fenhance_bins_meanCIstd(n_loc,1);
F0 = 1;
var0 = RC;
F1 = 2;
var1 = RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
xlabel('Flight Force','fontsize',10)
ylabel('pitch mid down-F mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% pitch mid upstroke
subplot(2,2,2)
n_loc = round((Rds_steady_meanCIstd(1)+1)*size(pitch_wb_L_bins,1)/2);
plot(F_mean_wb,pitch_wb_L_bins(n_loc,:),'.k')
plot(F_mean_wb,pitch_wb_R_bins(n_loc,:),'.k')
plot(F_Fenhance,pitch_wb_L_Fenhance_bins(n_loc,:),'.b')
plot(F_Fenhance,pitch_wb_R_Fenhance_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_Fenhance_bins_meanCIstd(n_loc,1)/Fenhance_norm;
F0 = 1;
var0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
F1 = 2;
var1 = var0 + RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
% xlabel('Flight Force','fontsize',10)
ylabel('pitch mid up','fontsize',10)

ang_min = 0;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,4)
plot(F_Fenhance,pitchMOD_wb_L_Fenhance_bins(n_loc,:),'.b')
plot(F_Fenhance,pitchMOD_wb_R_Fenhance_bins(n_loc,:),'.b')

% plot mean
RC = pitchMOD_wb_Fenhance_bins_meanCIstd(n_loc,1);
F0 = 1;
var0 = RC;
F1 = 2;
var1 = RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
xlabel('Flight Force','fontsize',10)
ylabel('pitch mid up-F mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 
