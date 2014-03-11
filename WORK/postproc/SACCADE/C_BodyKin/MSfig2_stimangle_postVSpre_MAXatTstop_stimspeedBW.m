subplot(3,3,6)
angle_pre_min = -180;
angle_pre_max = 0;
mean_min = -180;
mean_max = 0;

angle_pre = stim_angle_yaw_mirror_pre;
angle_post = stim_angle_vel_escape;
plotcolor = grey_color;
plot_angle_pre_post_circmean_extendedsection_meanlim_NOdp

angle_pre = stim_angle_yaw_mirror_pre(settings.expansion.speed==2);
angle_post = stim_angle_vel_escape(settings.expansion.speed==2);
plotcolor = [.8 .8 .8];
plot_angle_pre_post_circmean_extendedsection_meanlim_NOci

angle_pre = stim_angle_yaw_mirror_pre(settings.expansion.speed==1);
angle_post = stim_angle_vel_escape(settings.expansion.speed==1);
plotcolor = [.4 .4 .4];
plot_angle_pre_post_circmean_extendedsection_meanlim_NOci

angle_pre = stim_angle_yaw_mirror_pre(settings.expansion.speed==0);
angle_post = stim_angle_vel_escape(settings.expansion.speed==0);
plotcolor = 'w';
plot_angle_pre_post_circmean_extendedsection_meanlim_NOci

% plot([-180 0],[90 -90],'k--')
axis equal
xlabel('stimulus angle','fontsize',10) 
ylabel('escape heading','fontsize',10) 
set(gca,'xlim',[-180 0],'ylim',[-90 90])
set(gca,'XTick',[-180:90:180],'fontsize',8)
set(gca,'YTick',[-180:90:180],'fontsize',8) 
% grid on
