%% dev max

figure
subplot(3,2,1)
title('devmin@upstroke')
hold on
subplot(3,2,2)
title('devmin@upstroke mod coeff')
hold on
subplot(3,2,3)
hold on
subplot(3,2,4)
hold on
subplot(3,2,5)
hold on
subplot(3,2,6)
hold on

% min up dev bin location
start = 100;
stop = 200;
var_loc = dev_wb_steady_bins_meanCIstd(start:stop,1);
n_loc = start + find(var_loc == min(var_loc));

subplot(3,2,1)
% plot((ThirdMomentRatio),dev_wb_L_bins(n_loc,:),'.g')
plot((ThirdMomentRatio),dev_wb_R_bins(n_loc,:),'.g')
plot(ThirdMomentRatio_d3rdMom,dev_wb_R_d3rdMom_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_R_d3rdMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('ThirdMomentRatio','fontsize',10)
ylabel('right wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,3)
plot((ThirdMomentRatio),dev_wb_L_bins(n_loc,:),'.g')
% plot((ThirdMomentRatio),dev_wb_R_bins(n_loc,:),'.g')
plot((ThirdMomentRatio_d3rdMom),dev_wb_L_d3rdMom_bins(n_loc,:),'.b')

% plot trend
RC = devMOD_wb_L_d3rdMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = dev_wb_steady_bins_meanCIstd(n_loc,1);
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('ThirdMomentRatio','fontsize',10)
ylabel('left wing','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(3,2,5)
plot((ThirdMomentRatio),Ddev_wb_bins(n_loc,:),'.g')
% plot((ThirdMomentRatio),dev_wb_R_bins(n_loc,:),'.g')
plot((ThirdMomentRatio_d3rdMom),Ddev_wb_d3rdMom_bins(n_loc,:),'.b')

% plot trend
RC = DdevMOD_wb_d3rdMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = 0;
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('ThirdMomentRatio','fontsize',10)
ylabel('R - L','fontsize',10)

ang_min = -45;
ang_max = 45;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

%% coeff
% right
subplot(3,2,2)
plot((ThirdMomentRatio_d3rdMom),devMOD_wb_R_d3rdMom_bins(n_loc,:),'.b')

RC = devMOD_wb_R_d3rdMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('ThirdMomentRatio','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -180;
ang_max = 180;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% Left
subplot(3,2,4)
plot((ThirdMomentRatio_d3rdMom),devMOD_wb_L_d3rdMom_bins(n_loc,:),'.b')

RC = devMOD_wb_L_d3rdMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
% xlabel('ThirdMomentRatio','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -180;
ang_max = 180;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

% right-Left
subplot(3,2,6)
plot((ThirdMomentRatio_d3rdMom),DdevMOD_wb_d3rdMom_bins(n_loc,:),'.b')

RC = DdevMOD_wb_d3rdMom_bins_meanCIstd(n_loc,1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('ThirdMomentRatio','fontsize',10)
% ylabel('devmax mod coefficient','fontsize',10)

ang_min = -180;
ang_max = 180;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 


