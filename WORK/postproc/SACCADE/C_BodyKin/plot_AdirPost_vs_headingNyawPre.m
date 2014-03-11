subplot(1,2,1)
plotcolor = 'r'
angle_pre = heading_pre;
angle_post = Adir_post;
% plot_angle_pre_post_csaps
% plot_angle_pre_post_csaps_wrap
plot_angle_pre_post_circmean

xlabel('initial heading','fontsize',10) 
ylabel('escape A direction','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(1,2,2)
plotcolor = 'b'
angle_pre = yaw_pre;
angle_post = Adir_post;
% plot_angle_pre_post_csaps
% plot_angle_pre_post_csaps_wrap
plot_angle_pre_post_circmean

xlabel('initial yaw','fontsize',10) 
ylabel('escape A direction','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on
