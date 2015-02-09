
figure
subplot(3,1,1)
title('attitudes @ Amax')
hold on

% left turn
plotcolor = 'b'
angle_pre = stim_angle_yaw_mirror_pre;
angle_pre_min = -225;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;

subplot(3,1,1)
angle_post = roll_dot_dot_mirror_max;
plot_angle_pre_post_mean_extendedsection_meanlim
% xlabel('initial heading','fontsize',10) 
% axis equal
ylabel('roll accel','fontsize',10) 
set(gca,'ylim',[0 2e5])
set(gca,'YTick',[0:1e5:2e5],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on

subplot(3,1,2)
angle_post = pitch_dot_dot_max;
% plot_angle_pre_post_mean_mirrorsection
plot_angle_pre_post_mean_extendedsection_meanlim

% xlabel('initial heading','fontsize',10) 
% axis equal
ylabel('pitch accel','fontsize',10) 
set(gca,'ylim',[0 2e5])
set(gca,'YTick',[0:1e5:2e5],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal

subplot(3,1,3)
angle_post = yaw_dot_dot_mirror_max;
% plot_angle_pre_post_mean_doublemirrorsection
plot_angle_pre_post_mean_extendedsection_meanlim
% xlabel('initial heading','fontsize',10) 
% axis equal
ylabel('yaw accel','fontsize',10) 
set(gca,'ylim',[0 2e5])
set(gca,'YTick',[0:1e5:2e5],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal
