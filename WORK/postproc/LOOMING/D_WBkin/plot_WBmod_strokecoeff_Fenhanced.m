
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
plot(F_mean_wb,stroke_wb_L_bins(1,:),'.k')
plot(F_mean_wb,stroke_wb_R_bins(1,:),'.k')
plot(F_Fenhance,stroke_wb_L_Fenhance_bins(1,:),'.b')
plot(F_Fenhance,stroke_wb_R_Fenhance_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_Fenhance_bins_meanCIstd(1,1)/Fenhance_norm;
F0 = 1;
var0 = stroke_wb_steady_bins_meanCIstd(1,1);
F1 = 2;
var1 = var0 + RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
% xlabel('Flight Force','fontsize',10)
ylabel('stroke max','fontsize',10)

ang_min = 45;
ang_max = 135;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,3)
plot(F_Fenhance,strokeMOD_wb_L_Fenhance_bins(1,:),'.b')
plot(F_Fenhance,strokeMOD_wb_R_Fenhance_bins(1,:),'.b')

% plot mean
RC = strokeMOD_wb_Fenhance_bins_meanCIstd(1,1);
F0 = 1;
var0 = RC;
F1 = 2;
var1 = RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
xlabel('Flight Force','fontsize',10)
ylabel('strokemax-F mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% stroke min
n_loc = round(Rds_steady_meanCIstd(1)*size(pitch_wb_L_bins,1));

subplot(2,2,2)
plot(F_mean_wb,stroke_wb_L_bins(n_loc,:),'.k')
plot(F_mean_wb,stroke_wb_R_bins(n_loc,:),'.k')
plot(F_Fenhance,stroke_wb_L_Fenhance_bins(n_loc,:),'.b')
plot(F_Fenhance,stroke_wb_R_Fenhance_bins(n_loc,:),'.b')

% plot trend
RC = strokeMOD_wb_Fenhance_bins_meanCIstd(n_loc,1)/Fenhance_norm;
F0 = 1;
var0 = stroke_wb_steady_bins_meanCIstd(n_loc,1);
F1 = 2;
var1 = var0 + RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
% xlabel('Flight Force','fontsize',10)
ylabel('stroke min','fontsize',10)

ang_min = -90;
ang_max = 0;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,4)
plot(F_Fenhance,strokeMOD_wb_L_Fenhance_bins(n_loc,:),'.b')
plot(F_Fenhance,strokeMOD_wb_R_Fenhance_bins(n_loc,:),'.b')

% plot mean
RC = strokeMOD_wb_Fenhance_bins_meanCIstd(n_loc,1);
F0 = 1;
var0 = RC;
F1 = 2;
var1 = RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
xlabel('Flight Force','fontsize',10)
ylabel('strokemin-F mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 
