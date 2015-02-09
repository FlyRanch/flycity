figure
% subplot(3,1,2)
% title('mean attitudes')
% hold on

% left turn
plotcolor = 'b'
angle_pre = stim_angle_yaw_mirror_pre;
angle_pre_min = -225;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;


subplot(3,1,1)
angle_post = abs(dt_RollVel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('roll dt','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[0 .05])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[0 .05],'fontsize',8) 
grid on

subplot(3,1,2)
angle_post = abs(dt_YawVel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('yaw dt','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[0 .05])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[0 .05],'fontsize',8) 
grid on

subplot(3,1,3)
angle_post = abs(dt_PitchVel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('pitch dt','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[0 .05])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[0 .05],'fontsize',8) 
grid on


% right turn
plotcolor = 'r'
angle_pre = -stim_angle_yaw_mirror_pre;
angle_pre_min = -180;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

subplot(3,1,1)
angle_post = abs(dt_RollVel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,1,2)
angle_post = abs(dt_YawVel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,1,3)
angle_post = abs(dt_PitchVel_max_min);
plot_angle_pre_post_circmean_extendedsection_meanlim

