figure
angle_pre_min = -225;
angle_pre_max = 180;

subplot(2,2,1)
plotcolor = 'b'
angle_pre = stim_angle_vel_mirror_pre;
angle_post = stim_angle_accel_mirror_mean;
plot_angle_pre_post_circmean_extendedsection

xlabel('initial heading','fontsize',10) 
ylabel('A direction','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8) 
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(2,2,3)
plotcolor = 'b'
angle_pre = stim_angle_vel_mirror_pre;
angle_post = stim_angle_spn_mirror_mean;
plot_angle_pre_post_circmean_extendedsection

xlabel('initial heading','fontsize',10) 
ylabel('strokeplane normal','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(2,2,2)
plotcolor = 'b'
angle_pre = stim_angle_yaw_mirror_pre;
angle_post = stim_angle_accel_mirror_mean;
plot_angle_pre_post_circmean_extendedsection

xlabel('initial body angle','fontsize',10) 
ylabel('A direction','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(2,2,4)
plotcolor = 'b'
angle_pre = stim_angle_yaw_mirror_pre;
angle_post = stim_angle_spn_mirror_mean;
plot_angle_pre_post_circmean_extendedsection

xlabel('initial body angle','fontsize',10) 
ylabel('strokeplane normal','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

