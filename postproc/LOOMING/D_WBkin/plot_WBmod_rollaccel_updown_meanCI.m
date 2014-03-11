% plot wingbeat modifiactions for roll accel

t_bins = t_roll(:);
y_min = -45;
y_max = 45;

% stroke_up
var = Dstroke_rollaccel_up(:);
subplot(3,3,1)
hold on
plot_WBmeanCI

% pitch_up
var = Dpitch_rollaccel_up(:);
subplot(3,3,4)
hold on
plot_WBmeanCI

% dev_up
var = Ddev_rollaccel_up(:);
subplot(3,3,7)
hold on
plot_WBmeanCI

% stroke_down
var = Dstroke_rollaccel_down(:);
subplot(3,3,2)
hold on
plot_WBmeanCI

% pitch_down
var = Dpitch_rollaccel_down(:);
subplot(3,3,5)
hold on
plot_WBmeanCI


% dev_down
var = Ddev_rollaccel_down(:);
subplot(3,3,8)
hold on
plot_WBmeanCI

% dstroke
var = Ddstroke_rollaccel(:);
subplot(3,3,3)
hold on
plot_WBmeanCI

% dpitch
var = Ddpitch_rollaccel(:);
subplot(3,3,6)
hold on
plot_WBmeanCI


% ddev
var = Dddev_rollaccel(:);
subplot(3,3,9)
hold on
plot_WBmeanCI
