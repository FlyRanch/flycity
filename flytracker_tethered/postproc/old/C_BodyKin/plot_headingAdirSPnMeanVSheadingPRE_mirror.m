figure
angle_pre_min = -225;
angle_pre_max = 180;
plotcolor = 'b'

subplot(1,3,1)
angle_pre = stim_angle_vel_mirror_pre;
angle_post = stim_angle_vel_mirror_post;
plot_angle_pre_post_circmean_extendedsection

axis equal
xlabel('initial heading','fontsize',10) 
ylabel('escape heading','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180],'XTickLabel',[])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(1,3,2)
angle_pre = stim_angle_vel_mirror_pre;
angle_post = stim_angle_accel_mirror_mean;
plot_angle_pre_post_circmean_extendedsection

axis equal
xlabel('initial heading','fontsize',10) 
ylabel('A direction','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-180 180])
% set(gca,'XTick',[-180 -90 0 90 180],'XTickLabel',[])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(1,3,3)
angle_pre = stim_angle_vel_mirror_pre;
angle_post = stim_angle_spn_mirror_mean;
plot_angle_pre_post_circmean_extendedsection

axis equal
xlabel('initial heading','fontsize',10) 
ylabel('strokeplane normal','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on
