
figure
subplot(2,2,1)
hold on
subplot(2,2,2)
hold on
subplot(2,2,3)
hold on
subplot(2,2,4)
hold on

subplot(2,2,1)
var = f_wb_Fenhance;
h=hist(var);
hist(var)
axis([150 250 0 100*ceil(max(h)/100)])
    set(gca,'XTick',150:50:250) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('flap freq','fontsize',10)
%     grid on

subplot(2,2,2)
var = freqMOD_wb_Fenhance;
h=hist(var);
hist(var)
axis([-100 200 0 100*ceil(max(h)/100)])
    set(gca,'XTick',-200:100:400) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('f-F mod coefficient','fontsize',10)
%     grid on

subplot(2,2,3)
plot(F_mean_wb,f_wb_L,'.k')
plot(F_mean_wb,f_wb_R,'.k')
plot(F_Fenhance,f_wb_Fenhance,'.b')

% plot trend
RC = freqMOD_wb_Fenhance_meanCIstd(1)/Fenhance_norm;
F0 = 1;
f0 = f_wb_steady_meanCIstd(1);
F1 = 2;
f1 = f0 + RC;
plot([F0 F1],[f0 f1],'-r','linewidth',2)
xlabel('Flight Force','fontsize',10)
ylabel('flap freq','fontsize',10)

ang_min = 150;
ang_max = 250;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,4)
plot(F_Fenhance,freqMOD_wb_Fenhance,'.b')

% plot mean
RC = freqMOD_wb_Fenhance_meanCIstd(1);
F0 = 1;
f0 = RC;
F1 = 2;
f1 = RC;
plot([F0 F1],[f0 f1],'-r','linewidth',2)
xlabel('Flight Force','fontsize',10)
ylabel('f-F mod coefficient','fontsize',10)

ang_min = -100;
ang_max = 100;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

