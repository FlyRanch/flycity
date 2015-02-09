%% pitch mid down

figure
subplot(3,2,1)
title('pitch mid down')
hold on
subplot(3,2,2)
title('pitch middown mod coeff')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

% mid down bin location
n_loc = round(Rds_steady_meanCIstd(1) * size(pitch_wb_R_d2ndMom_bins,1)/2);

subplot(3,2,1)
% plot((SecondMomentRatio),pitch_wb_L_bins(n_loc,:),'.g')
plot((SecondMomentRatio),pitch_wb_R_bins(n_loc,:),'.g')
plot(SecondMomentRatio_d2ndMom,pitch_wb_R_d2ndMom_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_R_d2ndMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('SecondMomentRatio','fontsize',10)
ylabel('right wing','fontsize',10)

ang_min = 90;
ang_max = 180;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
plot((SecondMomentRatio),pitch_wb_L_bins(n_loc,:),'.g')
% plot((SecondMomentRatio),pitch_wb_R_bins(n_loc,:),'.g')
plot((SecondMomentRatio_d2ndMom),pitch_wb_L_d2ndMom_bins(n_loc,:),'.b')

% plot trend
RC = pitchMOD_wb_L_d2ndMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = pitch_wb_steady_bins_meanCIstd(n_loc,1);
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('SecondMomentRatio','fontsize',10)
ylabel('left wing','fontsize',10)

ang_min = 90;
ang_max = 180;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
plot((SecondMomentRatio),Dpitch_wb_bins(n_loc,:),'.g')
% plot((SecondMomentRatio),pitch_wb_R_bins(n_loc,:),'.g')
plot((SecondMomentRatio_d2ndMom),Dpitch_wb_d2ndMom_bins(n_loc,:),'.b')

% plot trend
RC = DpitchMOD_wb_d2ndMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = 0;
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('SecondMomentRatio','fontsize',10)
ylabel('R - L','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% coeff
% right
subplot(3,2,2)
plot((SecondMomentRatio_d2ndMom),pitchMOD_wb_R_d2ndMom_bins(n_loc,:),'.b')

RC = pitchMOD_wb_R_d2ndMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('SecondMomentRatio','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -90;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Left
subplot(3,2,4)
plot((SecondMomentRatio_d2ndMom),pitchMOD_wb_L_d2ndMom_bins(n_loc,:),'.b')

RC = pitchMOD_wb_L_d2ndMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('SecondMomentRatio','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -90;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% right-Left
subplot(3,2,6)
plot((SecondMomentRatio_d2ndMom),DpitchMOD_wb_d2ndMom_bins(n_loc,:),'.b')

RC = DpitchMOD_wb_d2ndMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('SecondMomentRatio','fontsize',10)
% ylabel('pitchmax mod coefficient','fontsize',10)

ang_min = -90;
ang_max = 90;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


