figure
% subplot(3,2,2)
% title('mean attitudes')
% hold on

% left turn
plotcolor = 'b'
angle_pre = stim_angle_yaw_mirror_pre;
angle_pre_min = -225;
angle_pre_max = 180;
mean_min = -180;
mean_max = 0;


subplot(3,2,1)
turn = turn_angle_vel_mirror;
% plot_angle_pre_turn_vectors_csaps
angle_post = turn;
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('heading turn','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-45 225])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(3,2,3)
angle_post = stim_angle_accel_mirror_mean;
angle_post = stim_angle_vel_mirror_post;
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('A direction','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(3,2,5)
angle_post = stim_angle_spn_mirror_mean;
plot_angle_pre_post_circmean_extendedsection_meanlim
ylabel('strokeplane normal','fontsize',10) 
set(gca,'xlim',[-225 45],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
set(gca,'YTick',[-180 -90 0 90 180],'fontsize',8) 
grid on

subplot(3,2,2)
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

subplot(3,2,4)
angle_post = calc_circ_mean_value(dpitch,n_pre,n_post);
plot_angle_pre_post_circmean_extendedsection_meanlim
% plot_angle_pre_post_circmean_mirrorsection
% xlabel('initial heading','fontsize',10) 
axis equal
ylabel('pitch','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal

subplot(3,2,6)
angle_post = calc_circ_mean_value(yaw_mirror,n_pre,n_post);
plot_angle_pre_post_circmean_extendedsection_meanlim
% plot_angle_pre_post_circmean_doublemirrorsection
% xlabel('initial heading','fontsize',10) 
axis equal
ylabel('yaw','fontsize',10) 
set(gca,'ylim',[-90 90])
set(gca,'YTick',[-180:45:180],'fontsize',8) 
set(gca,'xlim',[-225 45])
set(gca,'XTick',[-180:90:180],'fontsize',8) 
grid on
% axis equal

% right turn
plotcolor = 'r'
angle_pre = -stim_angle_yaw_mirror_pre;
angle_pre_min = -180;
angle_pre_max = 180;
mean_min = 0;
mean_max = 180;

subplot(3,2,1)
turn = -turn_angle_vel_mirror;
% plot_angle_pre_turn_vectors_csaps
angle_post = turn;
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'xlim',[-180 180],'ylim',[-180 180])
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
grid on

subplot(3,2,3)
angle_post = -stim_angle_accel_mirror_mean;
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)

subplot(3,2,5)
angle_post = -stim_angle_spn_mirror_mean;
plot_angle_pre_post_circmean_extendedsection_meanlim
set(gca,'XTick',[-180 -90 0 90 180],'fontsize',8)
xlabel('Stimulus angle','fontsize',10) 

subplot(3,2,2)
angle_post = -calc_circ_mean_value(roll_mirror,n_pre,n_post);
plot_angle_pre_post_circmean_extendedsection_meanlim

subplot(3,2,4)
angle_post = calc_circ_mean_value(dpitch,n_pre,n_post);
plot_angle_pre_post_circmean_extendedsection_meanlim
% set(gca,'xlim',[angle_pre_min angle_pre_max])

subplot(3,2,6)
angle_post = -calc_circ_mean_value(yaw_mirror,n_pre,n_post);
plot_angle_pre_post_circmean_extendedsection_meanlim
% set(gca,'xlim',[angle_pre_min angle_pre_max])
xlabel('Stimulus angle','fontsize',10) 
