% yaw rate VS DEV L-R
figure
subplot(3,2,1)
title('Downstroke L-R')
hold on
subplot(3,2,2)
title('Upstroke L-R')
hold on

subplot(3,2,1)
x = yaw_dot_mean_wb(:);
y = rad2deg(ddev_max_udsPREV(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
% xlabel('yaw rate')
ylabel('L-R @ max dev')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,2)
x = yaw_dot_mean_wb(:);
y = rad2deg(ddev_max_dus(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
% xlabel('yaw rate')
% ylabel('Ddev max')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,3)
x = yaw_dot_mean_wb(:);
y = rad2deg(ddev_min_ds(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
% xlabel('yaw rate')
ylabel('L-R @ min dev')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,4)
x = yaw_dot_mean_wb(:);
y = rad2deg(ddev_min_us(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
% xlabel('yaw rate')
% ylabel('Ddev min')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,5)
x = yaw_dot_mean_wb(:);
y = rad2deg(dAdev_ds(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('yaw rate')
ylabel('Adev diff')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,6)
x = yaw_dot_mean_wb(:);
y = rad2deg(dAdev_us(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('yaw rate')
% ylabel('dAdev')
axis([-2000 2500 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on
