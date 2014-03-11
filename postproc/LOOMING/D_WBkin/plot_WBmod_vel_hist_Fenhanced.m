
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
var = V_Fenhance;
h=hist(var);
hist(var)
axis([0 .7 0 100*ceil(max(h)/100)])
    set(gca,'XTick',150:50:250) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('speed','fontsize',10)
%     grid on

subplot(2,2,2)
var = velMOD_wb_Fenhance;
h=hist(var);
hist(var)
axis([-.5 1 0 100*ceil(max(h)/100)])
    set(gca,'XTick',-.5:.5:1) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('V-F mod coefficient','fontsize',10)
%     grid on

subplot(2,2,3)
plot(F_mean_wb,V_mean_wb,'.k')
plot(F_mean_wb,V_mean_wb,'.k')
plot(F_Fenhance,V_Fenhance,'.b')

% plot trend
RC = velMOD_wb_Fenhance_meanCIstd(1)/Fenhance_norm;
F0 = 1;
f0 = V_steady_meanCIstd(1);
F1 = 2;
f1 = f0 + RC;
plot([F0 F1],[f0 f1],'-r','linewidth',2)
xlabel('Flight Force','fontsize',10)
ylabel('speed','fontsize',10)

ang_min = 0;
ang_max = .7;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,4)
plot(F_Fenhance,velMOD_wb_Fenhance,'.b')

% plot mean
RC = velMOD_wb_Fenhance_meanCIstd(1);
F0 = 1;
f0 = RC;
F1 = 2;
f1 = RC;
plot([F0 F1],[f0 f1],'-r','linewidth',2)
xlabel('Flight Force','fontsize',10)
ylabel('f-F mod coefficient','fontsize',10)

ang_min = -.5;
ang_max = 1;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

