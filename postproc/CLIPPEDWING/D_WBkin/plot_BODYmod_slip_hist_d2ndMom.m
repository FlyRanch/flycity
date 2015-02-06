
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
ang_min = -180;
ang_max = 180;

subplot(2,2,1)
var = slip_d2ndMom;
h=hist(var);
hist(var)

axis([ang_min ang_max 0 100*ceil(max(h)/100)])
    set(gca,'XTick',ang_min:(ang_max-ang_min)/4:ang_max) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('slip','fontsize',10)
%     grid on

subplot(2,2,3)
plot(SecondMomentRatio,slip_mean_wb,'.g')
plot(SecondMomentRatio_d2ndMom,slip_d2ndMom,'.b')

% plot trend
RC = slipMOD_d2ndMom_meanCIstd(1);
x0 = xmax;
y0 = 0;
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('SecondMomentRatio','fontsize',10)
ylabel('slip','fontsize',10)

axis([xmin xmax ang_min ang_max])
set(gca,'XTick',xmin:(xmax-xmin)/2:xmax) 
set(gca,'YTick',ang_min:(ang_max-ang_min)/4:ang_max) 

%% gradient
ang_min = -2000;
ang_max = 2000;

subplot(2,2,2)
var = slipMOD_d2ndMom;
h=hist(var);
hist(var)

axis([ang_min ang_max 0 100*ceil(max(h)/100)])
    set(gca,'XTick',ang_min:(ang_max-ang_min)/4:ang_max) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('slip/d2ndMom (slipMOD)','fontsize',10)
%     grid on

subplot(2,2,4)
plot(SecondMomentRatio_d2ndMom,slipMOD_d2ndMom,'.b')

% plot mean
RC = slipMOD_d2ndMom_meanCIstd(1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('SecondMomentRatio','fontsize',10)
ylabel('slip/d2ndMom (slipMOD)','fontsize',10)

axis([xmin xmax ang_min ang_max])
set(gca,'XTick',xmin:(xmax-xmin)/2:xmax) 
set(gca,'YTick',ang_min:(ang_max-ang_min)/4:ang_max) 


