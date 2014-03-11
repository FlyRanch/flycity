%% stroke max

figure
subplot(3,2,1)
title('stroke max')
hold on
subplot(3,2,2)
title('strokemax mod coefficient')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

subplot(3,2,1)
% plot((M_R_mean_wb),stroke_wb_L_bins(1,:),'.k')
plot((M_R_mean_wb),stroke_wb_R_bins(1,:),'.k')
plot(M_R_TorqueAxisR,stroke_wb_R_TorqueAxisR_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_R_TorqueAxisR_bins_meanCIstd(1,1)/M_R_norm;
var0 = stroke_wb_steady_bins_meanCIstd(1,1);
x0 = 0;
x1 = max(M_R_TorqueAxisR);
var1 = var0 + RC * (x1-x0);
xmin1 = min(M_R_TorqueAxisR);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
ylabel('right wing','fontsize',10)

ang_min = 45;
ang_max = 135;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
plot((M_R_mean_wb),stroke_wb_L_bins(1,:),'.k')
% plot((M_R_mean_wb),stroke_wb_R_bins(1,:),'.k')
plot((M_R_TorqueAxisR),stroke_wb_L_TorqueAxisR_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_L_TorqueAxisR_bins_meanCIstd(1,1)/M_R_norm;
var0 = stroke_wb_steady_bins_meanCIstd(1,1);
x0 = 0;
x1 = max(M_R_TorqueAxisR);
var1 = var0 + RC * (x1-x0);
xmin1 = min(M_R_TorqueAxisR);
varmin1 = var0 + RC * (xmin1-x0);
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
ylabel('left wing','fontsize',10)

ang_min = 45;
ang_max = 135;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
plot((M_R_mean_wb),Dstroke_wb_bins(1,:),'.k')
% plot((M_R_mean_wb),stroke_wb_R_bins(1,:),'.k')
plot((M_R_TorqueAxisR),Dstroke_wb_TorqueAxisR_bins(1,:),'.b')

% plot trend
RC = DstrokeMOD_wb_TorqueAxisR_bins_meanCIstd(1,1)/M_R_norm;
var0 = 0;
x0 = 0;
x1 = max(M_R_TorqueAxisR);
var1 = var0 + RC * (x1-x0);
xmin1 = min(M_R_TorqueAxisR);
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
plot((M_R_TorqueAxisR),strokeMOD_wb_R_TorqueAxisR_bins(1,:),'.b')

RC = strokeMOD_wb_R_TorqueAxisR_bins_meanCIstd(1,1);
var0 = 0;
x0 = 0; 
x1 = max(M_R_TorqueAxisR); 
var1 = RC
xmin1 = min(M_R_TorqueAxisR); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Left
subplot(3,2,4)
plot((M_R_TorqueAxisR),strokeMOD_wb_L_TorqueAxisR_bins(1,:),'.b')

RC = strokeMOD_wb_L_TorqueAxisR_bins_meanCIstd(1,1);
var0 = 0;
x0 = 0; 
x1 = max(M_R_TorqueAxisR); 
var1 = RC
xmin1 = min(M_R_TorqueAxisR); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
% xlabel('rotR Torque','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% right-Left
subplot(3,2,6)
plot((M_R_TorqueAxisR),DstrokeMOD_wb_TorqueAxisR_bins(1,:),'.b')

RC = DstrokeMOD_wb_TorqueAxisR_bins_meanCIstd(1,1);
var0 = 0;
x0 = 0; 
x1 = max(M_R_TorqueAxisR); 
var1 = RC
xmin1 = min(M_R_TorqueAxisR); 
varmin1 = RC
plot([xmin1 x1],[varmin1 var1],'-r','linewidth',2)
xlabel('rotR Torque','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


