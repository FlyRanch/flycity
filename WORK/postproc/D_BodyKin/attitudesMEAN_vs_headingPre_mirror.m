figure
subplot(3,1,1)
title('mean attitudes')
hold on

% left turn
plotcolor = 'b'
angle_pre = stim_angle_vel_mirror_pre;
angle_pre_min = -225;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;

subplot(3,1,1)
angle_post = calc_circ_mean_value(roll_mirror,n_pre,n_post);
plot_angle_pre_post_circmean_extendedsection_meanlim
% xlabel('initial heading','fontsize',10) 
axis equal
ylabel('roll','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal

subplot(3,1,2)
angle_post = calc_circ_mean_value(pitch,n_pre,n_post);
plot_angle_pre_post_circmean_mirrorsection
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
angle_post = calc_circ_mean_value(yaw,n_pre,n_post);
plot_angle_pre_post_circmean_doublemirrorsection
% xlabel('initial heading','fontsize',10) 
axis equal
ylabel('yaw','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal