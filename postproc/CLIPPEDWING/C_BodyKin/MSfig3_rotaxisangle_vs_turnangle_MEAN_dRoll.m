% left turn

pitch_post = calc_circ_mean_value(dpitch,n_pre,n_post);
roll_post = calc_circ_mean_value(droll_mirror,n_pre,n_post);
% roll_post = calc_circ_mean_value(roll_mirror,n_pre,n_post);

% n_now = n_Amax;
% pitch_post = calc_value(dpitch,n_now);
% % roll_post = calc_value(droll_mirror,n_now);
% roll_post = calc_value(roll_mirror,n_now);

rot_dir = atand(pitch_post./roll_post)
angle_post = rot_dir;
% angle_post(angle_post<0) = abs(angle_post(angle_post<0))
% figure

subplot(3,3,6)
plotcolor = 'b'
plotcolor = [.5 .5 .5]
angle_pre = turn_angle_vel_mirror;
angle_pre_min = 0;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

plot_angle_pre_post_circmean_extendedsection_meanlim
% plot([-180 0],[90 -90],'k--')
% axis equal
xlabel('turn angle','fontsize',10) 
ylabel('mean rotation axis angle','fontsize',10) 
set(gca,'xlim',[angle_pre_min angle_pre_max],'ylim',[0 90])
set(gca,'XTick',[-180:45:180],'fontsize',8)
set(gca,'YTick',[-180:45:180],'fontsize',8) 
% grid on

