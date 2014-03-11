figure
angle_pre_min = -180;
angle_pre_max = 180;

subplot(2,2,1)
plotcolor = 'r'
angle_pre = stim_angle_vel_pre;
angle_post = stim_angle_accel_mean;
plot_angle_pre_post_circmean_doublemirrorsection

xlabel('heading pre','fontsize',10) 
ylabel('A direction','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8) 
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 

subplot(2,2,3)
plotcolor = 'b'
angle_pre = stim_angle_vel_pre;
angle_post = stim_angle_spn_mean;
plot_angle_pre_post_circmean_doublemirrorsection

xlabel('heading pre','fontsize',10) 
ylabel('body angle turn','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 

subplot(2,2,2)
plotcolor = 'r'
angle_pre = stim_angle_yaw_pre;
angle_post = stim_angle_accel_mean;
plot_angle_pre_post_circmean_doublemirrorsection

xlabel('initial body angle','fontsize',10) 
ylabel('heading turn','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 

subplot(2,2,4)
plotcolor = 'b'
angle_pre = stim_angle_yaw_pre;
angle_post = stim_angle_spn_mean;
plot_angle_pre_post_circmean_doublemirrorsection

xlabel('initial body angle','fontsize',10) 
ylabel('body angle turn','fontsize',10) 
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 

