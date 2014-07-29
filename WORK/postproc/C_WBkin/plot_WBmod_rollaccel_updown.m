% plot wingbeat modifiactions for roll accel
colmap = abs(colormap(gray)-1);

t_hist = t_roll(:);
t_bins = t_roll(:);
n_pol=9

biny_min = -180;
biny_max = 180;
biny = biny_min: (biny_max - biny_min)/ny :biny_max;

x_min = 0;
x_max = 1;

y_min = -45;
y_max = 45;

% stroke_up
var = Dstroke_rollaccel_up(:);

subplot(3,3,1)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing stroke','fontsize',10) 
title('roll accel mod left')

% pitch_up
var = Dpitch_rollaccel_up(:);

subplot(3,3,4)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing pitch','fontsize',10) 
% title('roll accel mod left')

% dev_up
var = Ddev_rollaccel_up(:);
subplot(3,3,7)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('deviation','fontsize',10) 
% title('roll accel left')

% stroke_down
var = Dstroke_rollaccel_down(:);
subplot(3,3,2)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing stroke','fontsize',10) 
title('roll accel mod right')

% pitch_down
var = Dpitch_rollaccel_down(:);
subplot(3,3,5)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing pitch','fontsize',10) 
% title('roll accel mod left')

% dev_down
var = Ddev_rollaccel_down(:);
subplot(3,3,8)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('deviation','fontsize',10) 
% title('roll accel left')


% dstroke
var = Ddstroke_rollaccel(:);
subplot(3,3,3)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing stroke','fontsize',10) 
title('roll accel mod right')

% dpitch
var = Ddpitch_rollaccel(:);
subplot(3,3,6)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('wing pitch','fontsize',10) 
% title('roll accel mod left')

% ddev
var = Dddev_rollaccel(:);
subplot(3,3,9)
plot_WBmean_heatmap
hold on
plot_WBmeanCI
% xlabel('time','fontsize',10) 
ylabel('deviation','fontsize',10) 
% title('roll accel left')


