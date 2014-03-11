% yaw rate VS PITCH L-R & STROKE L-R
figure
subplot(3,2,1)
title('wing stroke L-R')
hold on
subplot(3,2,2)
title('wing pitch L-R')
hold on

subplot(3,2,1)
x = yaw_dot_mean_wb(:);
y = rad2deg(dstroke_max_wb(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
% xlabel('yaw rate')
ylabel('L-R @ maximum')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,3)
x = yaw_dot_mean_wb(:);
y = rad2deg(dstroke_min_wb(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
% xlabel('yaw rate')
ylabel('L-R @ minimum')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,5)
x = yaw_dot_mean_wb(:);
y = rad2deg(dAstroke_wb(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('yaw rate')
ylabel('A diff')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,2)
x = yaw_dot_mean_wb(:);
y = rad2deg(dpitch_max_wb(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
% xlabel('yaw rate')
% ylabel('Dpitch max')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,4)
x = yaw_dot_mean_wb(:);
y = rad2deg(dpitch_min_wb(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
% xlabel('yaw rate')
% ylabel('Dpitch min')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,6)
x = yaw_dot_mean_wb(:);
y = rad2deg(dApitch_wb(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('yaw rate')
% ylabel('dApitch')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

