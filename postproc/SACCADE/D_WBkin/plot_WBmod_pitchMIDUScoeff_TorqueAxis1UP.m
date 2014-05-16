%% pitch mid up

figure
subplot(3,2,1)
title('pitch mid up')
hold on
subplot(3,2,2)
title('pitch midup mod coeff')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

% mid up bin location
n_loc = round((Rds_steady_meanCIstd(1)+1)*size(pitch_wb_R_TorqueAxis1UP_bins,1)/2);

subplot(3,2,1)
% plot((M_R_mean_wb),pitch_wb_L_bins(n_loc,:),'.k')
plot((M_R_mean_wb),pitch_wb_R_bins(n_loc,:),'.k')
plot(M_R_TorqueAxis1UP,pitch_wb_R_TorqueAxis1UP_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_R_TorqueAxis1UP_bins_meanCIstd(n_loc,1)/M_R_norm;
var0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0;
x1 = max(M_R_TorqueAxis1UP);
var1 = var0 + RC * (x1-x0);
xmin1 = min(M_R_TorqueAxis1UP);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
ylabel('right wing','fontsize',10)

ang_min = 0;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
plot((M_R_mean_wb),pitch_wb_L_bins(n_loc,:),'.k')
% plot((M_R_mean_wb),pitch_wb_R_bins(n_loc,:),'.k')
plot((M_R_TorqueAxis1UP),pitch_wb_L_TorqueAxis1UP_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_L_TorqueAxis1UP_bins_meanCIstd(n_loc,1)/M_R_norm;
var0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
x0 = 0;
x1 = max(M_R_TorqueAxis1UP);
var1 = var0 + RC * (x1-x0);
xmin1 = min(M_R_TorqueAxis1UP);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
ylabel('left wing','fontsize',10)

ang_min = 0;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
plot((M_R_mean_wb),Dpitch_wb_bins(n_loc,:),'.k')
% plot((M_R_mean_wb),pitch_wb_R_bins(n_loc,:),'.k')
plot((M_R_TorqueAxis1UP),Dpitch_wb_TorqueAxis1UP_bins(n_loc,:),'.b')

% plot trend
RC = DpitchMOD_wb_TorqueAxis1UP_bins_meanCIstd(n_loc,1)/M_R_norm;
var0 = 0;
x0 = 0;
x1 = max(M_R_TorqueAxis1UP);
var1 = var0 + RC * (x1-x0);
xmin1 = min(M_R_TorqueAxis1UP);
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
plot((M_R_TorqueAxis1UP),pitchMOD_wb_R_TorqueAxis1UP_bins(n_loc,:),'.b')

RC = pitchMOD_wb_R_TorqueAxis1UP_bins_meanCIstd(n_loc,1);
var0 = 0;
x0 = 0; 
x1 = max(M_R_TorqueAxis1UP); 
var1 = RC
xmin1 = min(M_R_TorqueAxis1UP); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Left
subplot(3,2,4)
plot((M_R_TorqueAxis1UP),pitchMOD_wb_L_TorqueAxis1UP_bins(n_loc,:),'.b')

RC = pitchMOD_wb_L_TorqueAxis1UP_bins_meanCIstd(n_loc,1);
var0 = 0;
x0 = 0; 
x1 = max(M_R_TorqueAxis1UP); 
var1 = RC
xmin1 = min(M_R_TorqueAxis1UP); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% right-Left
subplot(3,2,6)
plot((M_R_TorqueAxis1UP),DpitchMOD_wb_TorqueAxis1UP_bins(n_loc,:),'.b')

RC = DpitchMOD_wb_TorqueAxis1UP_bins_meanCIstd(n_loc,1);
var0 = 0;
x0 = 0; 
x1 = max(M_R_TorqueAxis1UP); 
var1 = RC
xmin1 = min(M_R_TorqueAxis1UP); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
xlabel('rotR Torque','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


