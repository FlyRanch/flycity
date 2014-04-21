angle_pre = -stim_angle_yaw_mirror_pre;
angle_pre_min = -180;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

subplot(3,1,1)
plotcolor = 'b'
angle_post = roll_dot_dot_mirror_max;
plot_angle_pre_post_mean_extendedsection_meanlim
plotcolor = 'r'
angle_post = roll_dot_dot_mirror_min;
plot_angle_pre_post_mean_extendedsection_meanlim

subplot(3,1,2)
plotcolor = 'b'
angle_post = pitch_dot_dot_max;
plot_angle_pre_post_mean_extendedsection_meanlim
plotcolor = 'r'
angle_post = pitch_dot_dot_min;
plot_angle_pre_post_mean_extendedsection_meanlim

subplot(3,1,3)
plotcolor = 'b'
angle_post = yaw_dot_dot_mirror_max;
plot_angle_pre_post_mean_extendedsection_meanlim
plotcolor = 'r'
angle_post = yaw_dot_dot_mirror_min;
plot_angle_pre_post_mean_extendedsection_meanlim

