
figure
subplot(2,2,1)
hold on
subplot(2,2,2)
hold on
subplot(2,2,3)
hold on
subplot(2,2,4)
hold on

%% max dev downstroke
var_loc = dev_wb_steady_bins_meanCIstd(1:50,1);
n_loc = find(var_loc == max(var_loc));

subplot(2,2,1)
plot(F_mean_wb,dev_wb_L_bins(n_loc,:),'.k')
plot(F_mean_wb,dev_wb_R_bins(n_loc,:),'.k')
plot(F_Fenhance,dev_wb_L_Fenhance_bins(n_loc,:),'.b')
plot(F_Fenhance,dev_wb_R_Fenhance_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_Fenhance_bins_meanCIstd(n_loc,1)/Fenhance_norm;
F0 = 1;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
F1 = 2;
var1 = var0 + RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
% xlabel('Flight Force','fontsize',10)
ylabel('max dev down','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,3)
plot(F_Fenhance,devMOD_wb_L_Fenhance_bins(n_loc,:),'.b')
plot(F_Fenhance,devMOD_wb_R_Fenhance_bins(n_loc,:),'.b')

% plot mean
RC = devMOD_wb_Fenhance_bins_meanCIstd(n_loc,1);
F0 = 1;
var0 = RC;
F1 = 2;
var1 = RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
xlabel('Flight Force','fontsize',10)
ylabel('max dev down mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% min dev mid
var_loc = dev_wb_steady_bins_meanCIstd(1:100,1);
n_loc = find(var_loc == min(var_loc));

subplot(2,2,2)
plot(F_mean_wb,dev_wb_L_bins(n_loc,:),'.k')
plot(F_mean_wb,dev_wb_R_bins(n_loc,:),'.k')
plot(F_Fenhance,dev_wb_L_Fenhance_bins(n_loc,:),'.b')
plot(F_Fenhance,dev_wb_R_Fenhance_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_Fenhance_bins_meanCIstd(n_loc,1)/Fenhance_norm;
F0 = 1;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
F1 = 2;
var1 = var0 + RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
% xlabel('Flight Force','fontsize',10)
ylabel('min dev down','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,4)
plot(F_Fenhance,devMOD_wb_L_Fenhance_bins(n_loc,:),'.b')
plot(F_Fenhance,devMOD_wb_R_Fenhance_bins(n_loc,:),'.b')

% plot mean
RC = devMOD_wb_Fenhance_bins_meanCIstd(n_loc,1);
F0 = 1;
var0 = RC;
F1 = 2;
var1 = RC;
plot([F0 F1],[var0 var1],'-r','linewidth',2)
xlabel('Flight Force','fontsize',10)
ylabel('min dev down-F mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 
