%% plot flightpath timeline tshift ON/OFF

close all
figure(1)
subplot(3,1,1)
plot(0,0,'r')
hold on
plot(0,0,'b')
legend('ON','OFF') 

for i=1:size(V,2)
    if settings.expansion.OFF(i) == 0
        color_now = 'r';
    elseif settings.expansion.OFF(i) == 1
        color_now = 'b';
    end

    subplot(3,2,1)
    hold on
    plot(t-t_shift(i),teta(:,i),color_now)
    subplot(3,2,3)
    hold on
    plot(t-t_shift(i),stim_angle_vel_plot(:,i),color_now)
    subplot(3,2,5)
    hold on
    plot(t-t_shift(i),stim_angle_yaw_plot(:,i),color_now)
    subplot(3,2,2)
    hold on
    plot(t-t_shift(i),V(:,i),color_now)
    subplot(3,2,4)
    hold on
    plot(t-t_shift(i),abs(An_hor(:,i)),color_now)
    subplot(3,2,6)
    hold on
    plot(t-t_shift(i),At_hor(:,i),color_now)
end



subplot(3,2,1)
axis([-nanmean(t_shift) .13-nanmean(t_shift) 0 180])
set(gca,'XTick',-.125:.025:.125)
set(gca,'YTick',0:90:180)
xlabel('time')
ylabel('optical angle')
grid on

subplot(3,2,3)
axis([-nanmean(t_shift) .13-nanmean(t_shift) -180 180])
set(gca,'XTick',-.125:.025:.125)
set(gca,'YTick',-180:180:180)
xlabel('time')
ylabel('heading angle')
grid on

subplot(3,2,5)
axis([-nanmean(t_shift) .13-nanmean(t_shift) -180 180])
set(gca,'XTick',-.125:.025:.125)
set(gca,'YTick',-180:180:180)
xlabel('time')
ylabel('yaw angle')
grid on

subplot(3,2,2)
axis([-nanmean(t_shift) .13-nanmean(t_shift) 0 1])
set(gca,'XTick',-.125:.025:.125)
set(gca,'YTick',0:.5:1)
xlabel('time')
ylabel('V')
grid on

subplot(3,2,4)
axis([-nanmean(t_shift) .13-nanmean(t_shift) 0 15])
set(gca,'XTick',-.125:.025:.125)
set(gca,'YTick',0:5:15)
xlabel('time')
ylabel('An')
grid on

subplot(3,2,6)
axis([-nanmean(t_shift) .13-nanmean(t_shift) -5 15])
set(gca,'XTick',-.125:.025:.125)
set(gca,'YTick',[-5;0;15])
xlabel('time')
ylabel('At')
grid on

