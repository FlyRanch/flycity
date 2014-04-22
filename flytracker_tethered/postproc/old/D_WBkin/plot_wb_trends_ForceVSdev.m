% normalized force VS DEV L&R
figure
subplot(3,2,1)
title('Downstroke Deviation')
hold on
subplot(3,2,2)
title('Upstroke Deviation')
hold on

subplot(3,2,1)
x = F_mean_wb(:);
y = rad2deg(dev_max_udsPREV_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(dev_max_udsPREV_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('normalized force')
ylabel('maximum dev')
axis([0.5 2.5 -45 45]) 
set(gca,'YTick',-90:45:180,'fontsize',8)
grid on

subplot(3,2,3)
x = F_mean_wb(:);
y = rad2deg(dev_min_ds_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(dev_min_ds_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('normalized force')
ylabel('minimum dev')
axis([0.5 2.5 -45 45]) 
set(gca,'YTick',-180:45:90,'fontsize',8)
grid on

subplot(3,2,5)
x = F_mean_wb(:);
y = rad2deg(Adev_ds_L(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(Adev_ds_R(:));
color_now = [.25 .25 .25];
plot_MAnCI_constwindow_circstd

xlabel('normalized force')
ylabel('Adev')
axis([0.5 2.5 -15 75]) 
set(gca,'YTick',-90:45:180,'fontsize',8)
grid on

subplot(3,2,2)
x = F_mean_wb(:);
y = rad2deg(dev_max_dus_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(dev_max_dus_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('normalized force')
% ylabel('Dpitch max')
axis([0.5 2.5 -45 45]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on

subplot(3,2,4)
x = F_mean_wb(:);
y = rad2deg(dev_min_us_L(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(dev_min_us_R(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd

% xlabel('normalized force')
% ylabel('Dpitch min')
axis([0.5 2.5 -45 45]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on

subplot(3,2,6)
x = F_mean_wb(:);
y = rad2deg(Adev_us_L(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd

x = F_mean_wb(:);
y = rad2deg(Adev_us_R(:));
color_now = [.25 .25 .25];
plot_MAnCI_constwindow_circstd

xlabel('normalized force')
% ylabel('dApitch')
axis([0.5 2.5 -15 75]) 
set(gca,'YTick',-180:45:180,'fontsize',8)
grid on

