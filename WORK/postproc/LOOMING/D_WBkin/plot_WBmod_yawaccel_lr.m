% plot wingbeat modifiactions for yaw accel
colmap = abs(colormap(gray)-1);

t_hist = t_yaw(:);
t_bins = t_yaw(:);

biny_min = -180;
biny_max = 180;
biny = biny_min: (biny_max - biny_min)/ny :biny_max;

x_min = 0;
x_max = 1;

y_min = -45;
y_max = 45;

% stroke_right
var = Dstroke_yawaccel_right(:);
subplot(3,3,1)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing stroke','fontsize',10) 
title('yaw accel mod left')

% pitch_right
var = Dpitch_yawaccel_right(:);
subplot(3,3,4)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing pitch','fontsize',10) 
% title('yaw accel mod left')

% dev_right
var = Ddev_yawaccel_right(:);
subplot(3,3,7)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('deviation','fontsize',10) 
% title('yaw accel left')

% stroke_left
var = Dstroke_yawaccel_left(:);
subplot(3,3,2)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing stroke','fontsize',10) 
title('yaw accel mod right')

% pitch_left
var = Dpitch_yawaccel_left(:);
subplot(3,3,5)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing pitch','fontsize',10) 
% title('yaw accel mod left')

% dev_left
var = Ddev_yawaccel_left(:);
subplot(3,3,8)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('deviation','fontsize',10) 
% title('yaw accel left')


% dstroke
var = Ddstroke_yawaccel(:);
subplot(3,3,3)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing stroke','fontsize',10) 
title('yaw accel mod right')

% dpitch
var = Ddpitch_yawaccel(:);
subplot(3,3,6)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing pitch','fontsize',10) 
% title('yaw accel mod left')

% ddev
var = Dddev_yawaccel(:);
subplot(3,3,9)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('deviation','fontsize',10) 
% title('yaw accel left')


