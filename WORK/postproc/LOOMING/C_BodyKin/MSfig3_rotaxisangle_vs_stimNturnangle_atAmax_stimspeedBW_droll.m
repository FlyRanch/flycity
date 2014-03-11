pitch_post = calc_circ_mean_value(dpitch,n_pre,n_post);
roll_post = calc_circ_mean_value(droll_mirror,n_pre,n_post);

n_now = n_Amax;
pitch_post = calc_value(dpitch,n_now);
roll_post = calc_value(droll_mirror,n_now);

% rot_dir = atand(pitch_post./roll_post)
rot_dir = atan2d(pitch_post,roll_post);

subplot(3,3,6)
angle_pre_min = -180;
angle_pre_max = 0;
mean_min = -180;
mean_max = 0;


angle_pre = stim_angle_yaw_mirror_pre;
angle_post = rot_dir;
plotcolor = grey_color;
plot_angle_pre_post_circmean_extendedsection_meanlim_NOdp

angle_pre = stim_angle_yaw_mirror_pre(settings.expansion.speed==2);
angle_post = rot_dir(settings.expansion.speed==2);
plotcolor = [.8 .8 .8];
plot_angle_pre_post_circmean_extendedsection_meanlim_NOci

angle_pre = stim_angle_yaw_mirror_pre(settings.expansion.speed==1);
angle_post = rot_dir(settings.expansion.speed==1);
plotcolor = [.4 .4 .4];
plot_angle_pre_post_circmean_extendedsection_meanlim_NOci

angle_pre = stim_angle_yaw_mirror_pre(settings.expansion.speed==0);
angle_post = rot_dir(settings.expansion.speed==0);
plotcolor = 'w';
plot_angle_pre_post_circmean_extendedsection_meanlim_NOci

plot([-180 0],[90 -90],'k--')
axis equal
xlabel('stimulus angle','fontsize',10) 
ylabel('rotation axis angle','fontsize',10) 
set(gca,'xlim',[-180 0],'ylim',[-45 90])
set(gca,'XTick',[-180:45:180],'fontsize',8)
set(gca,'YTick',[-180:45:180],'fontsize',8) 
