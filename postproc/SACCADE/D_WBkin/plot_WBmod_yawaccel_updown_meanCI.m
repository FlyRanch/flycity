% plot wingbeat modifiactions for yaw accel

t_bins = t_yaw(:);
y_min = -45;
y_max = 45;

% stroke_right
var = Dstroke_yawaccel_right(:);
subplot(3,3,1)
hold on
plot_WBmeanCI

% pitch_right
var = Dpitch_yawaccel_right(:);
subplot(3,3,4)
hold on
plot_WBmeanCI

% dev_right
var = Ddev_yawaccel_right(:);
subplot(3,3,7)
hold on
plot_WBmeanCI

% stroke_left
var = Dstroke_yawaccel_left(:);
subplot(3,3,2)
hold on
plot_WBmeanCI

% pitch_left
var = Dpitch_yawaccel_left(:);
subplot(3,3,5)
hold on
plot_WBmeanCI


% dev_left
var = Ddev_yawaccel_left(:);
subplot(3,3,8)
hold on
plot_WBmeanCI

% dstroke
var = Ddstroke_yawaccel(:);
subplot(3,3,3)
hold on
plot_WBmeanCI

% dpitch
var = Ddpitch_yawaccel(:);
subplot(3,3,6)
hold on
plot_WBmeanCI


% ddev
var = Dddev_yawaccel(:);
subplot(3,3,9)
hold on
plot_WBmeanCI
