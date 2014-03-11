angle_pre_min = -180;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

subplot(2,2,1)
plotcolor = 'r'
angle_pre = -stim_angle_vel_mirror_pre;
angle_post = -stim_angle_accel_mirror_mean;
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)

subplot(2,2,3)
plotcolor = 'r'
angle_pre = -stim_angle_vel_mirror_pre;
angle_post = -stim_angle_spn_mirror_mean;
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)

subplot(2,2,2)
plotcolor = 'r'
angle_pre = -stim_angle_yaw_mirror_pre;
angle_post = -stim_angle_accel_mirror_mean;
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)

subplot(2,2,4)
plotcolor = 'r'
angle_pre = -stim_angle_yaw_mirror_pre;
angle_post = -stim_angle_spn_mirror_mean;
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)

