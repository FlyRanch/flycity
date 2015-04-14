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
% plot((S2S3funcIntact),stroke_wb_L_bins(1,:),'.g')
plot((S2S3funcIntact),stroke_wb_R_bins(1,:),'.g')
plot(S2S3funcIntact_S2S3funcIntact,stroke_wb_R_S2S3funcIntact_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_R_S2S3funcIntact_bins_meanCIstd(1,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = stroke_wb_steady_bins_meanCIstd(1);
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3funcIntact','fontsize',10)
ylabel('Clipped wing','fontsize',10)

ang_min = 45;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
plot((S2S3funcIntact),stroke_wb_L_bins(1,:),'.g')
% plot((S2S3funcIntact),stroke_wb_R_bins(1,:),'.g')
plot((S2S3funcIntact_S2S3funcIntact),stroke_wb_L_S2S3funcIntact_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_L_S2S3funcIntact_bins_meanCIstd(1,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = stroke_wb_steady_bins_meanCIstd(1);
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3funcIntact','fontsize',10)
ylabel('Intact wing','fontsize',10)

ang_min = 45;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
plot((S2S3funcIntact),Dstroke_wb_bins(1,:),'.g')
% plot((S2S3funcIntact),stroke_wb_R_bins(1,:),'.g')
plot((S2S3funcIntact_S2S3funcIntact),Dstroke_wb_S2S3funcIntact_bins(1,:),'.b')

% plot trend
RC = DstrokeMOD_wb_S2S3funcIntact_bins_meanCIstd(1,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = 0;
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('S2S3funcIntact','fontsize',10)
ylabel('Clipped - Intact','fontsize',10)

ang_min = -15;
ang_max = 15;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% coeff
% Clipped
subplot(3,2,2)
plot((S2S3funcIntact_S2S3funcIntact),strokeMOD_wb_R_S2S3funcIntact_bins(1,:),'.b')

RC = strokeMOD_wb_R_S2S3funcIntact_bins_meanCIstd(1,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3funcIntact','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -1000;
ang_max = 1000;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Intact
subplot(3,2,4)
plot((S2S3funcIntact_S2S3funcIntact),strokeMOD_wb_L_S2S3funcIntact_bins(1,:),'.b')

RC = strokeMOD_wb_L_S2S3funcIntact_bins_meanCIstd(1,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3funcIntact','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -1000;
ang_max = 1000;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Clipped-Intact
subplot(3,2,6)
plot((S2S3funcIntact_S2S3funcIntact),DstrokeMOD_wb_S2S3funcIntact_bins(1,:),'.b')

RC = DstrokeMOD_wb_S2S3funcIntact_bins_meanCIstd(1,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('S2S3funcIntact','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -1000;
ang_max = 1000;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

