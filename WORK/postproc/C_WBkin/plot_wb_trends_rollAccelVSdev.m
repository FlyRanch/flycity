% ROLL acceleration VS DEV L-R
figure
subplot(3,2,1)
title('Downstroke L-R')
hold on
subplot(3,2,2)
title('Upstroke L-R')
hold on

% roll accel
x = roll_dot_dot_mean_wb(:);

subplot(3,2,1)
y = rad2deg(ddev_max_udsPREV(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
ylabel('max dev diff')
axis([-300000 150000 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on
fit_rollaccel_ddevmax_ds = MA;

subplot(3,2,2)
y = rad2deg(ddev_max_dus(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
% ylabel('Ddev max')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on
fit_rollaccel_ddevmax_us = MA;

subplot(3,2,3)
y = rad2deg(ddev_min_ds(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
ylabel('min dev diff')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on
fit_rollaccel_ddevmin_ds = MA;

subplot(3,2,4)
y = rad2deg(ddev_min_us(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
% xlabel('roll acceleration')
% ylabel('Ddev min')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on
fit_rollaccel_ddevmin_us = MA;

subplot(3,2,5)
y = rad2deg(dAdev_ds(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('roll acceleration')
ylabel('Adev diff')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on
fit_rollaccel_dAdev_ds = MA;

subplot(3,2,6)
y = rad2deg(dAdev_us(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('roll acceleration')
% ylabel('dAdev')
axis([-300000 150000 -45 45]) 
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on
fit_rollaccel_dAdev_us = MA;
