% right turn
plotcolor = 'r'
angle_pre = -stim_angle_vel_mirror_pre;
angle_pre_min = -180;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

subplot(3,1,1)
angle_post = -calc_circ_mean_value(roll_mirror,n_pre,n_post);
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,1,2)
set(gca,'xlim',[angle_pre_min angle_pre_max])
subplot(3,1,3)
set(gca,'xlim',[angle_pre_min angle_pre_max])
xlabel('initial heading','fontsize',10) 