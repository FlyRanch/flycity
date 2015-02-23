
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
ang_min = 150;
ang_max = 250;

subplot(2,2,1)
var = f_wb_d3rdMom;
h=hist(var);
hist(var)

axis([ang_min ang_max 0 100*ceil(max(h)/100)])
    set(gca,'XTick',ang_min:(ang_max-ang_min)/4:ang_max) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('freq','fontsize',10)
%     grid on

subplot(2,2,3)
plot(ThirdMomentRatio,f_wb_L,'.g')
plot(ThirdMomentRatio,f_wb_R,'.g')
plot(ThirdMomentRatio_d3rdMom,f_wb_d3rdMom,'.b')

% plot trend
RC = freqMOD_wb_d3rdMom_meanCIstd(1);
x0 = xmax;
y0 = f_wb_steady_meanCIstd(1);
x1 = xmin;
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('ThirdMomentRatio','fontsize',10)
ylabel('freq','fontsize',10)

axis([xmin xmax ang_min ang_max])
set(gca,'XTick',xmin:(xmax-xmin)/2:xmax) 
set(gca,'YTick',ang_min:(ang_max-ang_min)/4:ang_max) 

%% gradient
ang_min = -300;
ang_max = 300;

subplot(2,2,2)
var = freqMOD_wb_d3rdMom;
h=hist(var);
hist(var)

axis([ang_min ang_max 0 100*ceil(max(h)/100)])
    set(gca,'XTick',ang_min:(ang_max-ang_min)/4:ang_max) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('freq/d3rdMom (freqMOD)','fontsize',10)
%     grid on

subplot(2,2,4)
plot(ThirdMomentRatio_d3rdMom,freqMOD_wb_d3rdMom,'.b')

% plot mean
RC = freqMOD_wb_d3rdMom_meanCIstd(1);
x0 = xmax;
y0 = RC;
x1 = xmin;
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('ThirdMomentRatio','fontsize',10)
ylabel('freq/d3rdMom (freqMOD)','fontsize',10)

axis([xmin xmax ang_min ang_max])
set(gca,'XTick',xmin:(xmax-xmin)/2:xmax) 
set(gca,'YTick',ang_min:(ang_max-ang_min)/4:ang_max) 


