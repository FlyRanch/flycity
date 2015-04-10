%% dev max

figure
subplot(3,2,1)
title('devmax@upstroke')
hold on
subplot(3,2,2)
title('devmax@up mod coeff')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

start = 50;
stop = 150;
var_loc = dev_wb_steady_bins_meanCIstd(start:stop,1);
n_loc = start+find(var_loc == max(var_loc));

subplot(3,2,1)
% plot((S2S3AmpRatioFunc),dev_wb_L_bins(n_loc,:),'.g')
plot((S2S3AmpRatioFunc),dev_wb_R_bins(n_loc,:),'.g')
plot(S2S3AmpRatioFunc_S2S3AmpRatioFunc,dev_wb_R_S2S3AmpRatioFunc_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_R_S2S3AmpRatioFunc_bins_meanCIstd(n_loc,1);
x0 = S2S3AmpRatioFunc_NONclipped; 
y0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = xmax;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3AmpRatioFunc','fontsize',10)
ylabel('Clipped wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
plot((S2S3AmpRatioFunc),dev_wb_L_bins(n_loc,:),'.g')
% plot((S2S3AmpRatioFunc),dev_wb_R_bins(n_loc,:),'.g')
plot((S2S3AmpRatioFunc_S2S3AmpRatioFunc),dev_wb_L_S2S3AmpRatioFunc_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_L_S2S3AmpRatioFunc_bins_meanCIstd(n_loc,1);
x0 = S2S3AmpRatioFunc_NONclipped; 
y0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = xmax;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3AmpRatioFunc','fontsize',10)
ylabel('Intact wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
plot((S2S3AmpRatioFunc),Ddev_wb_bins(n_loc,:),'.g')
% plot((S2S3AmpRatioFunc),dev_wb_R_bins(n_loc,:),'.g')
plot((S2S3AmpRatioFunc_S2S3AmpRatioFunc),Ddev_wb_S2S3AmpRatioFunc_bins(n_loc,:),'.b')

% plot trend
RC = DdevMOD_wb_S2S3AmpRatioFunc_bins_meanCIstd(n_loc,1);
x0 = S2S3AmpRatioFunc_NONclipped; 
y0 = 0;
x1 = xmax;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('S2S3AmpRatioFunc','fontsize',10)
ylabel('Clipped - Intact','fontsize',10)

ang_min = -60;
ang_max = 30;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% coeff
% Clipped
subplot(3,2,2)
plot((S2S3AmpRatioFunc_S2S3AmpRatioFunc),devMOD_wb_R_S2S3AmpRatioFunc_bins(n_loc,:),'.b')

RC = devMOD_wb_R_S2S3AmpRatioFunc_bins_meanCIstd(n_loc,1);
x0 = S2S3AmpRatioFunc_NONclipped; 
y0 = RC;
x1 = xmax;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3AmpRatioFunc','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -0;
ang_max = 500;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Intact
subplot(3,2,4)
plot((S2S3AmpRatioFunc_S2S3AmpRatioFunc),devMOD_wb_L_S2S3AmpRatioFunc_bins(n_loc,:),'.b')

RC = devMOD_wb_L_S2S3AmpRatioFunc_bins_meanCIstd(n_loc,1);
x0 = S2S3AmpRatioFunc_NONclipped; 
y0 = RC;
x1 = xmax;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3AmpRatioFunc','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -500;
ang_max = 0;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Clipped-Intact
subplot(3,2,6)
plot((S2S3AmpRatioFunc_S2S3AmpRatioFunc),DdevMOD_wb_S2S3AmpRatioFunc_bins(n_loc,:),'.b')

RC = DdevMOD_wb_S2S3AmpRatioFunc_bins_meanCIstd(n_loc,1);
x0 = S2S3AmpRatioFunc_NONclipped; 
y0 = RC;
x1 = xmax;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('S2S3AmpRatioFunc','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -500;
ang_max = 0;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


