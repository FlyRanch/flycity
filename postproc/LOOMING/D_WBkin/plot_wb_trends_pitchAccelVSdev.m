% pitch acceleration VS DEV L&R
figure
subplot(3,2,1)
title('Downstroke Deviation')
hold on
subplot(3,2,2)
title('Upstroke Deviation')
hold on

% pitch accel
x = pitch_dot_dot_mean_wb(:);

subplot(3,2,1)
y = rad2deg(dev_max_udsPREV_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(dev_max_udsPREV_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

% xlabel('pitch acceleration')
ylabel('maximum dev')
axis([-200000 200000 -45 45]) 
set(gca,'YTick',-90:45:180,'fontsize',8)
grid on
fit_pitchaccel_devmax_ds = MA;

subplot(3,2,3)
y = rad2deg(dev_min_ds_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(dev_min_ds_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

% xlabel('pitch acceleration')
ylabel('minimum dev')
axis([-200000 200000 -45 45]) 
set(gca,'YTick',-180:45:90,'fontsize',8)
grid on
fit_pitchaccel_devmin_ds = MA;

subplot(3,2,5)
y = rad2deg(Adev_ds_L(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(Adev_ds_R(:));
color_now = [.25 .25 .25];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

xlabel('pitch acceleration')
ylabel('Adev')
axis([-200000 200000 -15 75]) 
set(gca,'YTick',-90:45:180,'fontsize',8)
grid on
fit_pitchaccel_Adev_ds = MA;

subplot(3,2,2)
y = rad2deg(dev_max_dus_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(dev_max_dus_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

% xlabel('pitch acceleration')
% ylabel('Dpitch max')
axis([-200000 200000 -45 45]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on
fit_pitchaccel_devmax_us = MA;

subplot(3,2,4)
y = rad2deg(dev_min_us_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(dev_min_us_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

% xlabel('pitch acceleration')
% ylabel('Dpitch min')
axis([-200000 200000 -45 45]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on
fit_pitchaccel_devmin_us = MA;

subplot(3,2,6)
y = rad2deg(Adev_us_L(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
MA_L = MA;

y = rad2deg(Adev_us_R(:));
color_now = [.25 .25 .25];
plot_MAnCI_constwindow_circstd
MA = mean([MA_L;MA]);

xlabel('pitch acceleration')
% ylabel('dApitch')
axis([-200000 200000 -15 75]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on
fit_pitchaccel_Adev_us = MA;

