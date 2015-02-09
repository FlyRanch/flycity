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
% plot((SecondMomentRatio),stroke_wb_L_bins(1,:),'.g')
plot((SecondMomentRatio),stroke_wb_R_bins(1,:),'.g')
plot(SecondMomentRatio_d2ndMom,stroke_wb_R_d2ndMom_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_R_d2ndMom_bins_meanCIstd(1,1);
x0 = xmax;
y0 = stroke_wb_steady_bins_meanCIstd(1);
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('SecondMomentRatio','fontsize',10)
ylabel('right wing','fontsize',10)

ang_min = 45;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
plot((SecondMomentRatio),stroke_wb_L_bins(1,:),'.g')
% plot((SecondMomentRatio),stroke_wb_R_bins(1,:),'.g')
plot((SecondMomentRatio_d2ndMom),stroke_wb_L_d2ndMom_bins(1,:),'.b')

% plot trend
RC = strokeMOD_wb_L_d2ndMom_bins_meanCIstd(1,1);
x0 = xmax;
y0 = stroke_wb_steady_bins_meanCIstd(1);
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('SecondMomentRatio','fontsize',10)
ylabel('left wing','fontsize',10)

ang_min = 45;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
plot((SecondMomentRatio),Dstroke_wb_bins(1,:),'.g')
% plot((SecondMomentRatio),stroke_wb_R_bins(1,:),'.g')
plot((SecondMomentRatio_d2ndMom),Dstroke_wb_d2ndMom_bins(1,:),'.b')

% plot trend
RC = DstrokeMOD_wb_d2ndMom_bins_meanCIstd(1,1);
x0 = xmax;
y0 = 0;
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('SecondMomentRatio','fontsize',10)
ylabel('R - L','fontsize',10)

ang_min = -15;
ang_max = 15;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% coeff
% right
subplot(3,2,2)
plot((SecondMomentRatio_d2ndMom),strokeMOD_wb_R_d2ndMom_bins(1,:),'.b')

RC = strokeMOD_wb_R_d2ndMom_bins_meanCIstd(1,1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('SecondMomentRatio','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -90;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Left
subplot(3,2,4)
plot((SecondMomentRatio_d2ndMom),strokeMOD_wb_L_d2ndMom_bins(1,:),'.b')

RC = strokeMOD_wb_L_d2ndMom_bins_meanCIstd(1,1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('SecondMomentRatio','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -90;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% right-Left
subplot(3,2,6)
plot((SecondMomentRatio_d2ndMom),DstrokeMOD_wb_d2ndMom_bins(1,:),'.b')

RC = DstrokeMOD_wb_d2ndMom_bins_meanCIstd(1,1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('SecondMomentRatio','fontsize',10)
% ylabel('strokemax mod coefficient','fontsize',10)

ang_min = -30;
ang_max = 30;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


