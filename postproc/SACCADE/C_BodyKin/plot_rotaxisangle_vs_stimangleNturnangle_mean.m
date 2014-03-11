% left turn

pitch_post = calc_circ_mean_value(dpitch,n_pre,n_post);
roll_post = calc_circ_mean_value(droll_mirror,n_pre,n_post);

% n_now = n_Amax;
% pitch_post = calc_value(dpitch,n_now);
% roll_post = calc_value(droll_mirror,n_now);

rot_dir = atand(pitch_post./roll_post)
% angle_post = abs(rot_dir);
angle_post = (rot_dir);

figure

subplot(1,2,1)
plotcolor = 'b'
angle_pre = stim_angle_yaw_mirror_pre;
angle_pre_min = -180;
angle_pre_max = 0;
mean_min = -180;
mean_max = 0;

plot_angle_pre_post_circmean_extendedsection_meanlim
axis equal
xlabel('stimulus angle','fontsize',10) 
ylabel('rotation axis angle','fontsize',10) 
set(gca,'xlim',[-180 0],'ylim',[-45 135])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on


subplot(1,2,2)
plotcolor = 'r'
angle_pre = turn_angle_vel_mirror;
angle_pre_min = 0;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

plot_angle_pre_post_circmean_extendedsection_meanlim
axis equal
xlabel('turn angle','fontsize',10) 
% ylabel('rotation axis angle','fontsize',10) 
set(gca,'xlim',[0 180],'ylim',[-45 135])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

