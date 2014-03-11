% yaw acceleration VS DEV L-R
figure
subplot(3,2,1)
title('Downstroke L-R')
hold on
subplot(3,2,2)
title('Upstroke L-R')
hold on

subplot(3,2,1)
plot(yaw_dot_dot_mean_wb(:),rad2deg(ddev_max_udsPREV(:)),'.')
% xlabel('yaw acceleration')
ylabel('Ddev max')
axis([-100000 200000 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)

subplot(3,2,2)
plot(yaw_dot_dot_mean_wb(:),rad2deg(ddev_max_dus(:)),'.')
% xlabel('yaw acceleration')
% ylabel('Ddev max')
axis([-100000 200000 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)

subplot(3,2,3)
plot(yaw_dot_dot_mean_wb(:),rad2deg(ddev_min_ds(:)),'.r')
% xlabel('yaw acceleration')
ylabel('Ddev min')
axis([-100000 200000 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)

subplot(3,2,4)
plot(yaw_dot_dot_mean_wb(:),rad2deg(ddev_min_us(:)),'.r')
% xlabel('yaw acceleration')
% ylabel('Ddev min')
axis([-100000 200000 -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)

subplot(3,2,5)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dAdev_ds(:)),'.k')
xlabel('yaw acceleration')
ylabel('dAdev')
axis([-100000 200000      -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)

subplot(3,2,6)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dAdev_us(:)),'.k')
xlabel('yaw acceleration')
% ylabel('dAdev')
axis([-100000 200000      -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)

% yaw acceleration VS PITCH L-R
figure
subplot(2,2,1)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dpitch_max_wb(:)),'.')
% xlabel('yaw acceleration')
ylabel('Dpitch max L-R')
axis([-100000 200000      -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)

subplot(2,2,2)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dpitch_min_wb(:)),'.')
% xlabel('yaw acceleration')
ylabel('Dpitch min L-R')
axis([-100000 200000      -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)

subplot(2,2,3)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dpitch_mean_wb(:)),'.k')
xlabel('yaw acceleration')
ylabel('dpitch mean L-R')
axis([-100000 200000      -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)

subplot(2,2,4)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dApitch_wb(:)),'.k')
xlabel('yaw acceleration')
ylabel('dApitch L-R')
axis([-100000 200000      -45 45])
set(gca,'YTick',-90:45:90,'fontsize',8)

% yaw VS STROKE L-R
figure
subplot(2,2,1)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dstroke_max_wb(:)),'.')
% xlabel('yaw acceleration')
ylabel('Dstroke max L-R')
axis([-100000 200000      -30 30])
set(gca,'YTick',-90:30:90,'fontsize',8)

subplot(2,2,2)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dstroke_min_wb(:)),'.')
% xlabel('yaw acceleration')
ylabel('Dstroke min L-R')
axis([-100000 200000      -30 30])
set(gca,'YTick',-90:30:90,'fontsize',8)

subplot(2,2,3)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dstroke_mean_wb(:)),'.k')
xlabel('yaw acceleration')
ylabel('dstroke mean L-R')
axis([-100000 200000      -30 30])
set(gca,'YTick',-90:30:90,'fontsize',8)

subplot(2,2,4)
plot(yaw_dot_dot_mean_wb(:),rad2deg(dAstroke_wb(:)),'.k')
xlabel('yaw acceleration')
ylabel('dAstroke L-R')
axis([-100000 200000      -30 30])
set(gca,'YTick',-90:30:90,'fontsize',8)

