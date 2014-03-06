n_now = n_Amax;

figure
subplot(3,1,1)
title('attitudes @ Amax')
hold on

% left turn
plotcolor = 'b'
angle_pre = stim_angle_vel_mirror_pre;
angle_pre_min = -225;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;

subplot(3,1,1)
angle_post = calc_value(roll_mirror,n_now);
plot_angle_pre_post_circmean_extendedsection_meanlim
% xlabel('initial heading','fontsize',10) 
axis equal
ylabel('roll','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on

subplot(3,1,2)
angle_post = calc_value(dpitch,n_now);
plot_angle_pre_post_circmean_extendedsection_meanlim
% xlabel('initial heading','fontsize',10) 
axis equal
ylabel('pitch','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal

subplot(3,1,3)
angle_post = calc_value(yaw_mirror,n_now);
plot_angle_pre_post_circmean_extendedsection_meanlim
% xlabel('initial heading','fontsize',10) 
axis equal
ylabel('yaw','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal
