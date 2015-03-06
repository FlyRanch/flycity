
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
ang_min = floor(10*min(vel_AeroCenterFuncs(:)))/10;
ang_max = ceil(10*max(vel_AeroCenterFuncs(:)))/10;
ang_min = 0;
ang_max = .5;

subplot(2,2,1)
var = vel_AeroCenterFuncs;
h=hist(var);
hist(var)

axis([ang_min ang_max 0 100*ceil(max(h)/100)])
    set(gca,'XTick',ang_min:(ang_max-ang_min)/4:ang_max) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('V','fontsize',10)
%     grid on

subplot(2,2,3)
plot(AeroCenterFuncRatio,V_mean_wb,'.g')
plot(AeroCenterFuncRatio_AeroCenterFuncs,vel_AeroCenterFuncs,'.b')

% plot trend
RC = velMOD_AeroCenterFuncs_meanCIstd(1);
x0 = AeroCenterFuncRatio_steady;
x1 = ceil(10*max(AeroCenterFuncRatio_AeroCenterFuncs(:)))/10;
y0 = V_steady_meanCIstd(1);
y1 = y0 + (x1-x0)*RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('AeroCenterFuncRatio','fontsize',10)
ylabel('V','fontsize',10)

axis([x0 x1 ang_min ang_max])
set(gca,'XTick',x0:(x1-x0)/2:x1) 
set(gca,'YTick',ang_min:(ang_max-ang_min)/4:ang_max) 

%% gradient
ang_min = -1;
ang_max = 1;

subplot(2,2,2)
var = velMOD_AeroCenterFuncs;
h=hist(var);
hist(var)

axis([ang_min ang_max 0 100*ceil(max(h)/100)])
    set(gca,'XTick',ang_min:(ang_max-ang_min)/4:ang_max) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('V/AeroCenterFuncRatio (velMOD)','fontsize',10)
%     grid on

subplot(2,2,4)
plot(AeroCenterFuncRatio_AeroCenterFuncs,velMOD_AeroCenterFuncs,'.b')

% plot mean
RC = velMOD_AeroCenterFuncs_meanCIstd(1);
plot([x0 x1],[RC RC],'-r','linewidth',2)
xlabel('AeroCenterFuncRatio','fontsize',10)
ylabel('V/AeroCenterFuncRatio (velMOD)','fontsize',10)

axis([x0 x1 ang_min ang_max])
set(gca,'XTick',x0:(x1-x0)/2:x1) 
set(gca,'YTick',ang_min:(ang_max-ang_min)/4:ang_max) 

