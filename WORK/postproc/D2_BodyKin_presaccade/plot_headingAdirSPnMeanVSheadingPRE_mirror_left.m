angle_pre_min = -180;
angle_pre_max = 180;
plotcolor = 'r'

subplot(1,3,1)
angle_pre = -stim_angle_vel_mirror_pre;
angle_post = -stim_angle_vel_mirror_post;
plot_angle_pre_post_circmean_extendedsection

subplot(1,3,2)
angle_pre = -stim_angle_vel_mirror_pre;
angle_post = -stim_angle_accel_mirror_mean;
plot_angle_pre_post_circmean_extendedsection

subplot(1,3,3)
angle_pre = -stim_angle_vel_mirror_pre;
angle_post = -stim_angle_spn_mirror_mean;
plot_angle_pre_post_circmean_extendedsection
