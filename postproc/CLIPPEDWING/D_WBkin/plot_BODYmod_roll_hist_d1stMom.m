
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
ang_min = -45;
ang_max = 45;

subplot(2,2,1)
var = roll_d1stMom;
h=hist(var);
hist(var)

axis([ang_min ang_max 0 100*ceil(max(h)/100)])
    set(gca,'XTick',ang_min:(ang_max-ang_min)/4:ang_max) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('roll','fontsize',10)
%     grid on

subplot(2,2,3)
plot(FirstMomentRatio,roll_mean_wb,'.g')
plot(FirstMomentRatio_d1stMom,roll_d1stMom,'.b')

% plot trend
RC = rollMOD_d1stMom_meanCIstd(1);
x0 = xmax;
y0 = 0;
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('FirstMomentRatio','fontsize',10)
ylabel('roll','fontsize',10)

axis([xmin xmax ang_min ang_max])
set(gca,'XTick',xmin:(xmax-xmin)/2:xmax) 
set(gca,'YTick',ang_min:(ang_max-ang_min)/4:ang_max) 

%% gradient
ang_min = -180;
ang_max = 180;

subplot(2,2,2)
var = rollMOD_d1stMom;
h=hist(var);
hist(var)

axis([ang_min ang_max 0 100*ceil(max(h)/100)])
    set(gca,'XTick',ang_min:(ang_max-ang_min)/4:ang_max) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('roll/d1stMom (rollMOD)','fontsize',10)
%     grid on

subplot(2,2,4)
plot(FirstMomentRatio_d1stMom,rollMOD_d1stMom,'.b')

% plot mean
RC = rollMOD_d1stMom_meanCIstd(1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('FirstMomentRatio','fontsize',10)
ylabel('roll/d1stMom (rollMOD)','fontsize',10)

axis([xmin xmax ang_min ang_max])
set(gca,'XTick',xmin:(xmax-xmin)/2:xmax) 
set(gca,'YTick',ang_min:(ang_max-ang_min)/4:ang_max) 


