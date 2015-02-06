angle_pre = stim_angle_yaw_mirror_pre;

t_min = t_stop - dt;
t_max = t_stop + dt;
t_window_subset = t_window(find(t_window>t_min & t_window<t_max));

angle_post = stim_angle_vel_plot(find(t_window>t_min & t_window<t_max));
% angle_post = stim_angle_vel_mirror_post;
% angle_post(angle_post<0) = abs(angle_post(angle_post<0))
% figure

subplot(3,3,9)
plotcolor = 'b'
plotcolor = [.5 .5 .5]
angle_pre_min = -180;
angle_pre_max = 0;
mean_min = -180;
mean_max = 0;

plot_angle_pre_post_circmean_extendedsection_meanlim
% plot([-180 0],[90 -90],'k--')
axis equal
xlabel('stimulus angle','fontsize',10) 
ylabel('escape heading','fontsize',10) 
set(gca,'xlim',[-180 0],'ylim',[-180 180])
set(gca,'XTick',[-180:90:180],'fontsize',8)
set(gca,'YTick',[-180:90:180],'fontsize',8) 
% grid on

