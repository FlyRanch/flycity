%% plot accel timeline tshift smooveOFF
figure
subplot(4,1,4)
plot(0,0,'r')
hold on
plot(0,0,'g')
plot(0,0,'b')
legend('165deg','64deg','64deg OFF') 

for i=1:size(u,2)
    if  settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 165
        subplot(4,1,1)
        hold on
        plot(t-t_shift(i),A(:,i),'r')
        subplot(4,1,2)
        hold on
        plot(t-t_shift(i),An(:,i),'r')
        subplot(4,1,3)
        hold on
        plot(t-t_shift(i),At(:,i),'r')
        subplot(4,1,4)
        hold on
        plot(t-t_shift(i),heading_dot(:,i),'r')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 64 && settings.expansion.OFF(i) == 0
        subplot(4,1,1)
        hold on
        plot(t-t_shift(i),A(:,i),'g')
        subplot(4,1,2)
        hold on
        plot(t-t_shift(i),An(:,i),'g')
        subplot(4,1,3)
        hold on
        plot(t-t_shift(i),At(:,i),'g')
        subplot(4,1,4)
        hold on
        plot(t-t_shift(i),heading_dot(:,i),'g')
    elseif settings.expansion.stepwise(i) == 0 && settings.expansion.speed(i) == 2 && settings.expansion.maxangle(i) == 64 && settings.expansion.OFF(i) == 1
        subplot(4,1,1)
        hold on
        plot(t-t_shift(i),A(:,i),'b')
        subplot(4,1,2)
        hold on
        plot(t-t_shift(i),An(:,i),'b')
        subplot(4,1,3)
        hold on
        plot(t-t_shift(i),At(:,i),'b')
        subplot(4,1,4)
        hold on
        plot(t-t_shift(i),heading_dot(:,i),'b')
    end
end

subplot(4,1,1)
% grid on
xlabel('time (s)')
ylabel('A')
% set(gca,'XTick',[-180 -90 0 90 180])
axis([t_start-nanmean(t_shift) t_stop-nanmean(t_shift) 0 30])

subplot(4,1,2)
% grid on
xlabel('time (s)')
ylabel('An')
axis([t_start-nanmean(t_shift) t_stop-nanmean(t_shift) 0 30])

subplot(4,1,3)
% grid on
xlabel('time (s)')
ylabel('At')
axis([t_start-nanmean(t_shift) t_stop-nanmean(t_shift) -10 20])

subplot(4,1,4)
% grid on
xlabel('time (s)')
ylabel('turn rate (deg/s)')
axis([t_start-nanmean(t_shift) t_stop-nanmean(t_shift) 0 10000])