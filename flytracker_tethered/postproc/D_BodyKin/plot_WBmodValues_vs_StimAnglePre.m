figure
% subplot(2,2,2)
% title('mean attitudes')
% hold on

angle_pre = stim_angle_yaw_mirror_pre;
angle_pre_min = -225;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;


subplot(2,2,1)
plotcolor = 'b'
angle_post = roll_dot_dot_norm_mirror_max;
plot_angle_pre_post_circmean_mirrorsection
plotcolor = 'r'
angle_post = roll_dot_dot_norm_mirror_min;
plot_angle_pre_post_circmean_mirrorsection
ylabel('MOD value roll accel','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-3 3])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-3 0 3],'fontsize',8) 
grid on

subplot(2,2,2)
plotcolor = 'b'
angle_post = pitch_dot_dot_norm_max;
plot_angle_pre_post_circmean_mirrorsection
plotcolor = 'r'
angle_post = pitch_dot_dot_norm_min;
plot_angle_pre_post_circmean_mirrorsection
ylabel('MOD value pitch accel','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-3 3])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-3 0 3],'fontsize',8) 
grid on

subplot(2,2,3)
plotcolor = 'b'
angle_post = yaw_dot_dot_norm_mirror_max;
plot_angle_pre_post_circmean_mirrorsection
plotcolor = 'r'
angle_post = yaw_dot_dot_norm_mirror_min;
plot_angle_pre_post_circmean_mirrorsection
ylabel('MOD value yaw accel','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-3 3])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-3 0 3],'fontsize',8) 
grid on

subplot(2,2,4)
plotcolor = 'b'
angle_post = F_norm_max;
plot_angle_pre_post_circmean_mirrorsection
ylabel('MOD value flight force','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-3 3])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-3 0 3],'fontsize',8) 
grid on










