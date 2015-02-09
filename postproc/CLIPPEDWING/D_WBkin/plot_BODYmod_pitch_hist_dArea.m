
figure
subplot(2,2,1)
hold on
subplot(2,2,2)
hold on
subplot(2,2,3)
hold on
subplot(2,2,4)
hold on

%% variable
ang_min = 15;
ang_max = 75;

subplot(2,2,1)
var = pitch_dArea;
h=hist(var);
hist(var)

axis([ang_min ang_max 0 100*ceil(max(h)/100)])
    set(gca,'XTick',ang_min:(ang_max-ang_min)/4:ang_max) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('pitch','fontsize',10)
%     grid on

subplot(2,2,3)
plot(AreaRatio,pitch_mean_wb,'.g')
plot(AreaRatio_dArea,pitch_dArea,'.b')

% plot trend
RC = pitchMOD_dArea_meanCIstd(1);
x0 = xmax;
y0 = pitch_body_steady;
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('AreaRatio','fontsize',10)
ylabel('pitch','fontsize',10)

axis([xmin xmax ang_min ang_max])
set(gca,'XTick',xmin:(xmax-xmin)/2:xmax) 
set(gca,'YTick',ang_min:(ang_max-ang_min)/4:ang_max) 

%% gradient
ang_min = -60;
ang_max = 60;

subplot(2,2,2)
var = pitchMOD_dArea;
h=hist(var);
hist(var)

axis([ang_min ang_max 0 100*ceil(max(h)/100)])
    set(gca,'XTick',ang_min:(ang_max-ang_min)/4:ang_max) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('pitch/dArea (pitchMOD)','fontsize',10)
%     grid on

subplot(2,2,4)
plot(AreaRatio_dArea,pitchMOD_dArea,'.b')

% plot mean
RC = pitchMOD_dArea_meanCIstd(1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('AreaRatio','fontsize',10)
ylabel('pitch/dArea (pitchMOD)','fontsize',10)

axis([xmin xmax ang_min ang_max])
set(gca,'XTick',xmin:(xmax-xmin)/2:xmax) 
set(gca,'YTick',ang_min:(ang_max-ang_min)/4:ang_max) 


