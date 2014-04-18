figure
subplot(3,3,1)
title('@ max roll accel')
hold on
subplot(3,3,2)
title('@ max yaw accel')
hold on
subplot(3,3,3)
title('@ max pitch accel')
hold on

%% left turn
angle_pre = stim_angle_yaw_mirror_pre;
angle_pre_min = -225;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;

% stroke
subplot(3,3,1)
plotcolor = 'b'
angle_post = Dstroke_RollAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Dstroke_RollAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('Dstroke','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

subplot(3,3,2)
plotcolor = 'b'
angle_post = Dstroke_YawAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Dstroke_YawAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim
% ylabel('Dstroke','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

subplot(3,3,3)
plotcolor = 'b'
angle_post = stroke_PitchAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = stroke_PitchAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim
% ylabel('Dstroke','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

% pitch
subplot(3,3,4)
plotcolor = 'b'
angle_post = Dpitch_RollAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Dpitch_RollAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('Dpitch','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

subplot(3,3,5)
plotcolor = 'b'
angle_post = Dpitch_YawAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Dpitch_YawAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim
% ylabel('Dpitch','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

subplot(3,3,6)
plotcolor = 'b'
angle_post = pitch_PitchAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = pitch_PitchAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim
% ylabel('Dpitch','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

% dev
subplot(3,3,7)
plotcolor = 'b'
angle_post = Ddev_RollAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Ddev_RollAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('Ddev','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

subplot(3,3,8)
plotcolor = 'b'
angle_post = Ddev_YawAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Ddev_YawAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim
% ylabel('Ddev','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

subplot(3,3,9)
plotcolor = 'b'
angle_post = dev_PitchAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = dev_PitchAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim
% ylabel('Ddev','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-30 30])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-30 0 30],'fontsize',8) 
grid on

%% right turn
angle_pre = -stim_angle_yaw_mirror_pre;
angle_pre_min = -180;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

% stroke
subplot(3,3,1)
plotcolor = 'b'
angle_post = Dstroke_RollAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Dstroke_RollAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,3,2)
plotcolor = 'b'
angle_post = Dstroke_YawAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Dstroke_YawAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,3,3)
plotcolor = 'b'
angle_post = stroke_PitchAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = stroke_PitchAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim

% pitch
subplot(3,3,4)
plotcolor = 'b'
angle_post = Dpitch_RollAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Dpitch_RollAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,3,5)
plotcolor = 'b'
angle_post = Dpitch_YawAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Dpitch_YawAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,3,6)
plotcolor = 'b'
angle_post = pitch_PitchAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = pitch_PitchAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim

% dev
subplot(3,3,7)
plotcolor = 'b'
angle_post = Ddev_RollAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Ddev_RollAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,3,8)
plotcolor = 'b'
angle_post = Ddev_YawAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = Ddev_YawAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,3,9)
plotcolor = 'b'
angle_post = dev_PitchAccel_max;
plot_angle_pre_post_circmean_extendedsection_meanlim
plotcolor = 'r'
angle_post = dev_PitchAccel_min;
plot_angle_pre_post_circmean_extendedsection_meanlim

