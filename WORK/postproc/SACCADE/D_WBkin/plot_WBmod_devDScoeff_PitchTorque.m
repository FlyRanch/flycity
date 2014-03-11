
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
plot(Mpitch_mean_wb,dev_wb_L_bins(n_loc,:),'.k')
plot(Mpitch_mean_wb,dev_wb_R_bins(n_loc,:),'.k')
plot(Mpitch_PitchTorque,dev_wb_L_PitchTorque_bins(n_loc,:),'.b')
plot(Mpitch_PitchTorque,dev_wb_R_PitchTorque_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_PitchTorque_bins_meanCIstd(n_loc,1)/Mpitch_norm;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0;
x1 = max(Mpitch_PitchTorque);
var1 = var0 + RC * (x1-x0);
xmin1 = min(Mpitch_PitchTorque);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('Pitch Torque','fontsize',10)
ylabel('max dev down','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,3)
plot(Mpitch_PitchTorque,devMOD_wb_L_PitchTorque_bins(n_loc,:),'.b')
plot(Mpitch_PitchTorque,devMOD_wb_R_PitchTorque_bins(n_loc,:),'.b')

% plot mean
RC = devMOD_wb_PitchTorque_bins_meanCIstd(n_loc,1);
x0 = min(Mpitch_PitchTorque); 
x1 = max(Mpitch_PitchTorque); 
var0 = RC;
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Pitch Torque','fontsize',10)
ylabel('max dev down-F mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% min dev mid
var_loc = dev_wb_steady_bins_meanCIstd(1:100,1);
n_loc = find(var_loc == min(var_loc));

subplot(2,2,2)
plot(Mpitch_mean_wb,dev_wb_L_bins(n_loc,:),'.k')
plot(Mpitch_mean_wb,dev_wb_R_bins(n_loc,:),'.k')
plot(Mpitch_PitchTorque,dev_wb_L_PitchTorque_bins(n_loc,:),'.b')
plot(Mpitch_PitchTorque,dev_wb_R_PitchTorque_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_PitchTorque_bins_meanCIstd(n_loc,1)/Mpitch_norm;
var0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0;
x1 = max(Mpitch_PitchTorque);
var1 = var0 + RC * (x1-x0);
xmin1 = min(Mpitch_PitchTorque);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('Pitch Torque','fontsize',10)
ylabel('min dev down','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,4)
plot(Mpitch_PitchTorque,devMOD_wb_L_PitchTorque_bins(n_loc,:),'.b')
plot(Mpitch_PitchTorque,devMOD_wb_R_PitchTorque_bins(n_loc,:),'.b')

% plot mean
RC = devMOD_wb_PitchTorque_bins_meanCIstd(n_loc,1);
x0 = min(Mpitch_PitchTorque); 
x1 = max(Mpitch_PitchTorque); 
var0 = RC;
var1 = RC;
plot([x0 x1],[var0 var1],'-r','linewidth',2)
xlabel('Pitch Torque','fontsize',10)
ylabel('min dev down-F mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 
