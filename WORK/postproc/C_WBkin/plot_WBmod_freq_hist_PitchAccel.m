
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
var = f_wb_PitchAccel;
h=hist(var);
hist(var)
axis([150 250 0 100*ceil(max(h)/100)])
    set(gca,'XTick',150:50:250) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('flap freq','fontsize',10)
%     grid on

subplot(2,2,2)
var = freqMOD_wb_PitchAccel;
h=hist(var);
hist(var)
axis([-200 200 0 100*ceil(max(h)/100)])
    set(gca,'XTick',-1e4:1e4:1e4) 
%     set(gca,'YTick',0:100:500,'fontsize',8)
    xlabel('f mod coefficient','fontsize',10)
%     grid on

subplot(2,2,3)
plot(pitch_dot_dot_mean_wb,f_wb_L,'.k')
plot(pitch_dot_dot_mean_wb,f_wb_R,'.k')
plot(pitchaccel_PitchAccel,f_wb_PitchAccel,'.b')

% plot trend
RC = freqMOD_wb_PitchAccel_meanCIstd(1)/pitchaccel_norm;
x0 = 0;
y0 = f_wb_steady_meanCIstd(1);
x1 = max(pitchaccel_PitchAccel);
y1 = y0 + RC * (x1-x0);
xmin1 = min(pitchaccel_PitchAccel);
ymin1 = y0 + RC * (xmin1-x0);
plot([x0 x1],[y0 y1],'-r','linewidth',2)
plot([x0 xmin1],[y0 ymin1],'-r','linewidth',2)
xlabel('Pitch Accel','fontsize',10)
ylabel('flap freq','fontsize',10)

ang_min = 150;
ang_max = 250;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 

subplot(2,2,4)
plot(pitchaccel_PitchAccel,freqMOD_wb_PitchAccel,'.b')

% plot mean
RC = freqMOD_wb_PitchAccel_meanCIstd(1);
x0 = min(pitchaccel_PitchAccel);
y0 = RC;
x1 = max(pitchaccel_PitchAccel);
y1 = RC;
plot([x0 x1],[y0 y1],'-r','linewidth',2)
xlabel('Pitch Accel','fontsize',10)
ylabel('f mod coefficient','fontsize',10)


ang_min = -100;
ang_max = 100;
axis([xmin xmax ang_min ang_max])
set(gca,'YTick',ang_min:(ang_max-ang_min)/2:ang_max) 
