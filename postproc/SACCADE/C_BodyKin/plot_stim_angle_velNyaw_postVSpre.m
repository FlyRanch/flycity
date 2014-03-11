subplot(2,2,1)
plotcolor = 'r'
angle_pre = stim_angle_vel_pre;
angle_post = stim_angle_vel_post;
% plot_angle_pre_post_csaps
% plot_angle_pre_post_csaps_wrap
plot_angle_pre_post_circmean

% xlabel('initial heading','fontsize',10) 
ylabel('escape heading','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180]) 
set(gca,'XTick',[-180 -90 0 90 180],'xticklabel',[])
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(2,2,3)
plotcolor = 'b'
angle_pre = stim_angle_vel_pre;
angle_post = stim_angle_yaw_post;
% plot_angle_pre_post_csaps
% plot_angle_pre_post_csaps_wrap
plot_angle_pre_post_circmean

xlabel('initial heading','fontsize',10) 
ylabel('escape yaw','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180]) 
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(2,2,2)
plotcolor = 'r'
angle_pre = stim_angle_yaw_pre;
angle_post = stim_angle_vel_post;
% plot_angle_pre_post_csaps
% plot_angle_pre_post_csaps_wrap
plot_angle_pre_post_circmean

% xlabel('initial yaw','fontsize',10) 
% ylabel('escape heading','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180]) 
set(gca,'XTick',[-180 -90 0 90 180],'xticklabel',[])
set(gca,'YTick',[-180 -90 0 90 180],'yticklabel',[],'fontsize',8) 
grid on

subplot(2,2,4)
plotcolor = 'b'
angle_pre = stim_angle_yaw_pre;
angle_post = stim_angle_yaw_post;
% plot_angle_pre_post_csaps
% plot_angle_pre_post_csaps_wrap
plot_angle_pre_post_circmean

xlabel('initial yaw','fontsize',10) 
% ylabel('escape yaw','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180]) 
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180],'yticklabel',[],'fontsize',8) 
grid on

