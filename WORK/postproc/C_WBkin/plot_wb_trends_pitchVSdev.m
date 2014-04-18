% pitch VS DEV L&R
figure
subplot(3,2,1)
title('Downstroke Deviation')
hold on
subplot(3,2,2)
title('Upstroke Deviation')
hold on

subplot(3,2,1)
x = pitch_mean_wb(:);
y = rad2deg(dev_max_udsPREV_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = pitch_mean_wb(:);
y = rad2deg(dev_max_udsPREV_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('pitch')
ylabel('maximum angle')
axis([-15 105 -45 45]) 
set(gca,'XTick',-90:45:180,'fontsize',8)
set(gca,'YTick',-90:45:180,'fontsize',8)
grid on

subplot(3,2,3)
x = pitch_mean_wb(:);
y = rad2deg(dev_min_ds_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = pitch_mean_wb(:);
y = rad2deg(dev_min_ds_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('pitch')
ylabel('minimum angle')
axis([-15 105 -45 45]) 
set(gca,'XTick',-90:45:180,'fontsize',8)
set(gca,'YTick',-180:45:90,'fontsize',8)
grid on

subplot(3,2,5)
x = pitch_mean_wb(:);
y = rad2deg(Adev_ds_L(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd

x = pitch_mean_wb(:);
y = rad2deg(Adev_ds_R(:));
color_now = [.25 .25 .25];
plot_MAnCI_constwindow_circstd

xlabel('pitch')
ylabel('Amplitude')
axis([-15 105 -15 75]) 
set(gca,'XTick',-90:45:180,'fontsize',8)
set(gca,'YTick',-90:45:180,'fontsize',8)
grid on

subplot(3,2,2)
x = pitch_mean_wb(:);
y = rad2deg(dev_max_dus_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = pitch_mean_wb(:);
y = rad2deg(dev_max_dus_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('pitch')
% ylabel('Dpitch max')
axis([-15 105 -45 45]) 
set(gca,'XTick',-90:45:180,'fontsize',8)
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on

subplot(3,2,4)
x = pitch_mean_wb(:);
y = rad2deg(dev_min_us_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = pitch_mean_wb(:);
y = rad2deg(dev_min_us_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('pitch')
% ylabel('Dpitch min')
axis([-15 105 -45 45]) 
set(gca,'XTick',-90:45:180,'fontsize',8)
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on

subplot(3,2,6)
x = pitch_mean_wb(:);
y = rad2deg(Adev_us_L(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd

x = pitch_mean_wb(:);
y = rad2deg(Adev_us_R(:));
color_now = [.25 .25 .25];
plot_MAnCI_constwindow_circstd

xlabel('pitch')
% ylabel('dApitch')
axis([-15 105 -15 75]) 
set(gca,'XTick',-90:45:180,'fontsize',8)
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on

