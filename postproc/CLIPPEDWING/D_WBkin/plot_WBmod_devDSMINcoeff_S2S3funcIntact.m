%% dev max

figure
subplot(3,2,1)
title('devmin@downstroke')
hold on
subplot(3,2,2)
title('devmin@down mod coeff')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

% min down dev bin location
var_loc = dev_wb_steady_bins_meanCIstd(1:100,1);
n_loc = find(var_loc == min(var_loc));

subplot(3,2,1)
% plot((S2S3funcIntact),dev_wb_L_bins(n_loc,:),'.g')
plot((S2S3funcIntact),dev_wb_R_bins(n_loc,:),'.g')
plot(S2S3funcIntact_S2S3funcIntact,dev_wb_R_S2S3funcIntact_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_R_S2S3funcIntact_bins_meanCIstd(n_loc,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3funcIntact','fontsize',10)
ylabel('Clipped wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
plot((S2S3funcIntact),dev_wb_L_bins(n_loc,:),'.g')
% plot((S2S3funcIntact),dev_wb_R_bins(n_loc,:),'.g')
plot((S2S3funcIntact_S2S3funcIntact),dev_wb_L_S2S3funcIntact_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_L_S2S3funcIntact_bins_meanCIstd(n_loc,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3funcIntact','fontsize',10)
ylabel('Intact wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
plot((S2S3funcIntact),Ddev_wb_bins(n_loc,:),'.g')
% plot((S2S3funcIntact),dev_wb_R_bins(n_loc,:),'.g')
plot((S2S3funcIntact_S2S3funcIntact),Ddev_wb_S2S3funcIntact_bins(n_loc,:),'.b')

% plot trend
RC = DdevMOD_wb_S2S3funcIntact_bins_meanCIstd(n_loc,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = 0;
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('S2S3funcIntact','fontsize',10)
ylabel('Clipped - Intact','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% coeff
% Clipped
subplot(3,2,2)
plot((S2S3funcIntact_S2S3funcIntact),devMOD_wb_R_S2S3funcIntact_bins(n_loc,:),'.b')

RC = devMOD_wb_R_S2S3funcIntact_bins_meanCIstd(n_loc,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3funcIntact','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -2000;
ang_max = 2000;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Intact
subplot(3,2,4)
plot((S2S3funcIntact_S2S3funcIntact),devMOD_wb_L_S2S3funcIntact_bins(n_loc,:),'.b')

RC = devMOD_wb_L_S2S3funcIntact_bins_meanCIstd(n_loc,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('S2S3funcIntact','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -2000;
ang_max = 2000;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Clipped-Intact
subplot(3,2,6)
plot((S2S3funcIntact_S2S3funcIntact),DdevMOD_wb_S2S3funcIntact_bins(n_loc,:),'.b')

RC = DdevMOD_wb_S2S3funcIntact_bins_meanCIstd(n_loc,1);
x0 = S2S3funcIntact_NONclipped; 
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('S2S3funcIntact','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = 0;
ang_max = 4000;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


