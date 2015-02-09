figure
dm=20
dn=20
plotcolor = 'b'
n_now = n_Amax;
angle_pre = stim_angle_yaw_pre;
angle_pre_min = -180;
angle_pre_max = 180;

subplot(3,1,1)
angle_post = calc_value(roll,n_now);
plot_angle_pre_post_circmean_doublemirrorsection
% xlabel('initial heading','fontsize',10) 
ylabel('roll','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on
% axis equal
title('attitudes @ Amax')

subplot(3,1,2)
angle_post = calc_value(pitch,n_now);
plot_angle_pre_post_circmean_mirrorsection
% xlabel('initial heading','fontsize',10) 
ylabel('pitch','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on
% axis equal

subplot(3,1,3)
angle_post = calc_value(yaw,n_now);
plot_angle_pre_post_circmean_doublemirrorsection
% xlabel('initial heading','fontsize',10) 
ylabel('yaw','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
grid on
% axis equal
