figure
angle_pre_min = -225;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;

subplot(2,2,1)
plotcolor = 'b'
angle_pre = stim_angle_vel_mirror_pre;
turn = turn_angle_vel_mirror;
% plot_angle_pre_turn_vectors_csaps
angle_post = turn;
plot_angle_pre_post_circmean_extendedsection_meanlim

xlabel('heading pre','fontsize',10) 
ylabel('heading turn','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-45 225])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8) 
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(2,2,3)
plotcolor = 'b'
angle_pre = stim_angle_vel_mirror_pre;
turn = turn_angle_yaw_mirror
% plot_angle_pre_turn_vectors_csaps
angle_post = turn;
plot_angle_pre_post_circmean_extendedsection_meanlim

xlabel('heading pre','fontsize',10) 
ylabel('body angle turn','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-45 225])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(2,2,2)
plotcolor = 'b'
angle_pre = stim_angle_yaw_mirror_pre;
turn = turn_angle_vel_mirror;
% plot_angle_pre_turn_vectors_csaps
angle_post = turn;
plot_angle_pre_post_circmean_extendedsection_meanlim

xlabel('initial body angle','fontsize',10) 
ylabel('heading turn','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-45 225])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(2,2,4)
plotcolor = 'b'
angle_pre = stim_angle_yaw_mirror_pre;
turn = turn_angle_yaw_mirror
% plot_angle_pre_turn_vectors_csaps
angle_post = turn;
plot_angle_pre_post_circmean_extendedsection_meanlim

xlabel('initial body angle','fontsize',10) 
ylabel('body angle turn','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-45 225])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

