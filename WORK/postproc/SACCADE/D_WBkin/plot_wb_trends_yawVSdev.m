% yaw VS DEV L-R
figure
subplot(3,2,1)
title('Downstroke L-R')
hold on
subplot(3,2,2)
title('Upstroke L-R')
hold on

subplot(3,2,1)
x = yaw_mean_wb(:);
y = rad2deg(ddev_max_udsPREV(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
% xlabel('yaw')
ylabel('L-R @ max dev')
axis([-45 90 -45 45])
set(gca,'XTick',-90:45:90,'fontsize',8)
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,2)
x = yaw_mean_wb(:);
y = rad2deg(ddev_max_dus(:));
color_now = [1 0 0];
plot_MAnCI_constwindow_circstd
% xlabel('yaw')
% ylabel('Ddev max')
axis([-45 90 -45 45]) 
set(gca,'XTick',-90:45:90,'fontsize',8)
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,3)
x = yaw_mean_wb(:);
y = rad2deg(ddev_min_ds(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
% xlabel('yaw')
ylabel('L-R @ min dev')
axis([-45 90 -45 45]) 
set(gca,'XTick',-90:45:90,'fontsize',8)
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,4)
x = yaw_mean_wb(:);
y = rad2deg(ddev_min_us(:));
color_now = [0 0 1];
plot_MAnCI_constwindow_circstd
plot(xp,yp,'color',.5*color_now,'linewidth',3)
% xlabel('yaw')
% ylabel('Ddev min')
axis([-45 90 -45 45]) 
set(gca,'XTick',-90:45:90,'fontsize',8)
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,5)
x = yaw_mean_wb(:);
y = rad2deg(dAdev_ds(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('yaw')
ylabel('Adev diff')
axis([-45 90 -45 45]) 
set(gca,'XTick',-90:45:90,'fontsize',8)
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on

subplot(3,2,6)
x = yaw_mean_wb(:);
y = rad2deg(dAdev_us(:));
color_now = [.5 .5 .5];
plot_MAnCI_constwindow_circstd
xlabel('yaw')
% ylabel('dAdev')
axis([-45 90 -45 45]) 
set(gca,'XTick',-90:45:90,'fontsize',8)
set(gca,'YTick',-90:45:90,'fontsize',8)
grid on